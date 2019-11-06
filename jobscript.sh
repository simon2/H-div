#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10034b
#QSUB -ug gr10034
#QSUB -W 1:00
#QSUB -A p=1:t=36:c=36:m=120G
#============ Shell Script ===========

mpiexec.hydra ./hmat_div_cilk 36
