#### Information ####
# Conducts PC-RELATE for GENESIS analysis based on PLINK Binaries
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-05-26

#### Command line input ####
library(optparse)
optlist = list(
  make_option(c('-i', '--in'), dest = 'bfile', help = 'PLINK binary prefix'),
  make_option(c('-f','--force'), dest = 'force', action = 'store_true', default = F,
              help = 'force overwrite')
)
args = parse_args(OptionParser(option_list = optlist))
print('Input options')
print(args)

#### Main execution block ####
main = function(args){
  library(SNPRelate)
  library(GWASTools)
  library(GENESIS)
  
  #### Convert input PLINK binary to gds format ####
  gdsfile = paste0(args$bfile,'.gds')
  if (! file.exists(gdsfile) | args$force) snpgdsBED2GDS(
    bed.fn = paste0(args$bfile,'.bed'), bim.fn = paste0(args$bfile, '.bim'),
    fam.fn = paste0(args$bfile,'.fam'), out.gdsfn = gdsfile
  )
  
  #### LD Pruning ####
  gds = snpgdsOpen(gdsfile)
  snpfile = paste0(args$bfile,'.snpset.rdata')
  if (! file.exists(snpfile) | args$force){
    snpset = snpgdsLDpruning(gds, method = 'corr', slide.max.bp = 1e6,
                           ld.threshold = sqrt(.1), verbose = T, num.thread = 8)
    save(snpset, file = snpfile)
  } else load(snpfile)
  pruned = unlist(snpset, use.names = F)
  
  #### Kinship Matrix ####
  kinfile = paste0(args$bfile,'.king.rdata')
  if (! file.exists(kinfile) | args$force){
    kin = snpgdsIBDKING(gds, num.thread = 8)
    eid = kin$sample.id
    kinmat = kin$kinship; colnames(kinmat) = eid; rownames(kinmat) = eid
    save(kinmat, file = kinfile)
  } else load(kinfile)
  snpgdsClose(gds)
  
  #### PC-AiR ####
  pcairfile = paste0(args$bfile,'.pcair.rdata')
  geno = GdsGenotypeReader(gdsfile)
  genodata = GenotypeData(geno)
  if (! file.exists(pcairfile) | args$force){
    pcairres = pcair(genodata, kinobj = kinmat, divobj = kinmat, 
                     snp.include = pruned, num.cores = 8)
    save(pcairres, file = pcairfile)
  } else load(pcairfile)
  pdf(paste0(args$bfile,'.pcair.pdf'))
  plot(pcairres); dev.off()
  
  #### PC-RELATE ####
  genoiter = GenotypeBlockIterator(genodata, snpInclude = pruned)
  pcrelatefile = paste0(args$bfile,'.pcrelate.rdata')
  if (! file.exists(pcrelatefile) | args$force){
    pcrelateres = pcrelate(genoiter, pcs = pcairres$vectors[,1:2], 
                           sample.block.size = 10000,
                           BPPARAM = BiocParallel::MulticoreParam(workers = 4))
    pcrelategrm = pcrelateToMatrix(pcrelateres)
    save(pcrelateres, pcrelategrm, file = pcrelatefile)
  }
}
main(args)
warnings()