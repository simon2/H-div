#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10034b
#QSUB -ug gr10034
#QSUB -W 4:00
#QSUB -A p=1:t=1:c=1:m=120G
#============ Shell Script ============

mpiexec.hydra ./hd
