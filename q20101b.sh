#!/bin/bash
#============ PBS Options ============
#QSUB -q gr20100b
#QSUB -ug gr20115
#QSUB -W 4:00
#QSUB -A p=1:t=36:c=36:m=120G
#============ Shell Script ===========

mpiexec.hydra ./hdc 36
