CC=icc
CXX=icpc
COPT=-std=c99 -O3 -xavx2 -g
COPT_OMP=-qopenmp

TARGETS = \
hmat_div \
hmat_div_direct \
hmat_div_cilk \
hmat_div_omp \
hmat_div_BCT_cilk_list_reducer \
hmat_div_BCT_cilk_lock \
hmat_div_BCT_cilk_malloc \
hmat_div_CT_cilk_parByLevel \
hmat_div.cpp \
hmat.c \
hmat_div_tcell \
hmat_div_locality \
hmat_dist

.phony: all

all: $(TARGETS)

# Serial version
hmat_div: hmat_div.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^

# Serial version with direct operations on data but not on data index
hmat_div_direct: hmat_div_direct.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^

# Cilk Plus version
hmat_div_cilk: hmat_div_cilk.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^

# Parallelized by OpenMP (only for experiments)
hmat_div_omp: hmat_div_omp.c
	$(CC) $(COPT) $(COPT_OMP) -o $@ $<

# Parallelized BCT by list_reducer in Cilk Plus for C++
hmat_div_BCT_cilk_list_reducer: hmat_div_BCT_cilk_list_reducer.cpp
	$(CXX) $(COPT) -o $@ $<

# Parallelized BCT by Cilk Plus and use mutex lock for result recording.
hmat_div_BCT_cilk_lock: hmat_div_BCT_cilk_lock.c
	$(CC) $(COPT) -o $@ $<

# Parallelized BCT by Cilk Plus and do malloc if origin space is not enough for result recording.
hmat_div_BCT_cilk_malloc: hmat_div_block_cilk_malloc.c
	$(CC) $(COPT) -o $@ $<

# Parallelized CT by Cilk Plus and control parallization by tree level. 
hmat_div_CT_cilk_parByLevel: hmat_div_CT_cilk_parByLevel.c
	$(CC) $(COPT) -o $@ $<

# Sequential C++ version
hmat_div.cpp: hmat_div.cpp
	$(CXX) $(COPT) -o $@ $<

# Sequential (SC version)
hmat.c: hmat.sc
	sc2c $<

# Parallelized by Tascell using same methods as hmat_div_cilk.c
hmat_div_tcell: hmat_div.tcell
	sc2c $<

# Parallelized by Tascell and taking data locality into consider.
hmat_div_locality: hmat_div_locality.tcell
	sc2c $<

# Parallelized by Tascell on distributed memory systems (baseline version)
hmat_dist: hmat_dist.tcell
	sc2c $<
