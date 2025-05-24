#!/bin/bash
source /home/yh464/.bashrc
cd /rds/project/rds-Nl99R8pHODQ/ref
conda activate wd

python variant_ref.py ${SLURM_ARRAY_TASK_ID}
