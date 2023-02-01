CC=icc
CXX=icpc
CMPI=mpiicc
COPT=-std=c11 -mcmodel=medium -shared-intel -O3 -xavx2 -g
COPT_OMP=-qopenmp
COPT_MKL=-mkl
COPT_MKL_SEQ=-mkl=sequential

TARGETS = \
hmat_div \
hmat_div_direct \
hmat_div_array \
hmat_div_cilk \
hmat_div_omp \
hmat_div_BCT_cilk_list_reducer \
hmat_div_BCT_cilk_lock \
hmat_div_BCT_cilk_malloc \
hmat_div_CT_cilk_parByLevel \
hmat_div_cpp \
hmat_array_filling \
hmat_array_filling_wBCT \
hmat_array_filling_MPI \
hmat_array_filling_dynamic \
hmat_sc

.phony: all

all: $(TARGETS)

# Serial version
hmat_div: hmat_div.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^

# Serial version with direct operations on data but not on data index
hmat_div_direct: hmat_div_direct.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^

hmat_div_array: hmat_div_array.c data/bem_file.c
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
hmat_div_cpp: hmat_div.cpp
	$(CXX) $(COPT) -o $@ $<

# Filling versions
hmat_array_filling: hmat_array_filling.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^ $(COPT_MKL_SEQ)

hmat_array_filling_wBCT: hmat_array_filling_wBCT.c data/bem_file.c
	$(CC) $(COPT) -o $@ $^ $(COPT_MKL_SEQ)

# Filling in parallel by MPI & OpenMP
hmat_array_filling_MPI: hmat_array_filling_MPI.c data/bem_file.c
	$(CMPI) $(COPT) -o $@ $^ $(COPT_MKL)

hmat_array_filling_dynamic: hmat_array_filling_dynamic.c data/bem_file.c
	$(CMPI) $(COPT) $(COPT_OMP) -o $@ $^ $(COPT_MKL_SEQ)

# Sequential (SC version)
hmat_sc: hmat.sc
	sc2c $<

# Parallelized by Tascell using same methods as hmat_div_cilk.c
hmat_div_tcell: hmat_div.tcell
	sc2c $<

# Parallelized by Tascell and taking data locality into consider.
hmat_div_locality: hmat_div_locality.tcell
	sc2c $<
