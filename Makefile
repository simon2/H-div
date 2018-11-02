CC=icc
COPT=-std=c99 -O3 -xavx2
COPT_OMP=-qopenmp

TARGETS = \
hacapk_division \
hacapk_division_cilk \
hacapk_division_cilk_lock \
hacapk_division_cilk_merge \
hacapk_division_cilk_merge_v1 \
hacapk_division_cilk_merge_v2 \
hacapk_division_omp

.phony: all

all: $(TARGETS)

# Changed the nesting order of the parallel loops for finding min/max.
# ("loop over particles -> loop over dimensions" in this version)
hacapk_division_cilk_merge: hacapk_division_cilk_merge.c
	$(CC) $(COPT) -o $@ $<

# In each step of constructing a cluster tree, find min/max and sort particles in parallel.
hacapk_division_cilk_merge_v2: hacapk_division_cilk_merge_v2.c
	$(CC) $(COPT) -o $@ $<

# Used merging strategy to resolve contention 
hacapk_division_cilk_merge_v1: hacapk_division_cilk_merge_v1.c
	$(CC) $(COPT) -o $@ $<

# Replaced three spawn calls to (nested) cilk_for loops
hacapk_division_cilk: hacapk_division_cilk.c
	$(CC) $(COPT) -o $@ $<

# Using pthread_mutex_lock to reslove contention when writing leaf elements of BC tree.
hacapk_division_cilk_lock: hacapk_division_cilk_lock.c
	$(CC) $(COPT) -o $@ $<

# Parallelized the 1st level of cnst. BC tree by OpenMP (only for experiments)
hacapk_division_omp: hacapk_division_omp.c
	$(CC) $(COPT) $(COPT_OMP) -o $@ $<

# Sequential
hacapk_division: hacapk_division.c
	$(CC) $(COPT) -o $@ $<

# Sequential (SC version)
hmat.c: hmat.sc
	sc2c $<
hmat: hmat.c
	$(CC) $(COPT) -o $@ $<