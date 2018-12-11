#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10034b
#QSUB -ug gr10034
#QSUB -W 0:10
#QSUB -A p=1:t=36:c=36:m=120G
#============ Shell Script ===========

mpiexec.hydra ./hdcm 36
