# Synopsis

H-div is a parallelized implementation of construction of hierarchical matrices using Cilk Plus and Tascell based on a sequential Fortran implementation in HACApK library coded by Prof. Akihiro Ida.

- To learn more about HACApK, please visit [HACApK introduction](http://ppopenhpc.cc.u-tokyo.ac.jp/ppopenhpc/2017/01/31/ppopen-applbem-ver-0-5-0/) or [HACApK github repository](https://github.com/Post-Peta-Crest/ppOpenHPC/tree/MATH/HACApK).

- To learn more about Tascell, please visit [Tascell introduction](http://ais.sys.i.kyoto-u.ac.jp/~task/tascell/index.html) or [Tascell github repository](https://github.com/tascell/sc-tascell).

- For details of this implementation on shared memory systems, please read paper [_Parallelization of Matrix Partitioning in Construction of Hierarchical Matrices using Task Parallel Languages_](https://www.jstage.jst.go.jp/article/ipsjjip/27/0/27_840/_article/-char/ja/).

- For details of this implementation on distributed memory systems, please read paper [_Parallelization of Matrix Partitioning in Hierarchical Matrix Construction on Distributed Memory Systems_](https://www.jstage.jst.go.jp/article/ipsjjip/30/0/30_742/_article/-char/ja/).

- For details of parallel implementation with both matrix partitioning and filling operation using Tascell, please read paper [Construction of Hierarchical Matrix on Distributed Memory Systems using a Task Parallel Language](https://ieeexplore.ieee.org/abstract/document/10062615).

# Requirements
- An Intel multi-core CPU
- To run Cilk Plus versions, Intel C++ Compiler version >= 17

  (Note: Cilk Plus is no longer in support from Intel and will be removed from Intel compiler some day in the future.)
- To run tascell versions on shared memory systems
  - tascell compiler version later than Jan 21, 2019.
  - GCC version >= 4.8.5 (or ICC with compatibility of GCC version higher than 4.8.5)
- To run tascell versions on distributed memory systems
  - tascell compiler version later than May 15, 2022 of branch ``mpi-bcst``.
  - Intel C++ Compiler version >= 17
  - Intel MPI version >= 17

# Installation
```
git clone https://github.com/simon2/H-div.git
```

# Compile & Execution
```
make hmat_div
./hmat_div
```

# File Explanation
## 1. matrix partitioning only
1. Sequential
    - hmat_div.c: The first sequential implementation using C directly translated from Fortran implementation in HACApK.
    - hmat_div_direct.c: Sequential C implementation but exchange data elements directly instead of exchange their index. This is the baseline of sequential implementation in paper.
    - hmat_div.cpp: Sequential implementation of C++.
    - hmat_div.sc: Sequential implementation using S-expression-based syntax. This is the base of Tascell.
    - hmat_div_array.c: CT is not in linked-tree manner, but use a pre-allocated array. (will be discribed in next paper)
2. Cilk Plus
    - hmat_div_cilk.c: Final implementation using Cilk Plus.
    - hmat_div_BCT_cilk_list_reducer.cpp: Tried CILK_LIST_REDUCER, but not included in paper due to bad performance.
    - hmat_div_CT_cilk_parByLevel.c: Tried to swith parallel code to sequential code by tree level information. Not included in paper due to bad performance.
    - hmat_div_BCT_cilk_malloc.c: Create private BCT array for each worker by using malloc function. Not included in paper due to bad performance.
3. Tascell
    - hmat_div.tcell: The final implementation using Tascell.
    - hmat_div_locality.tcell: Tried to make upper tree levels execute sequentially and execute in parallel in lower levels, for better data locality. However, the result of speedup is not essential.
    - hmat_dist.tcell: Baseline version parallelized on distributed memory systems.
    - hmat_dist_bcst.tcell: Add broad-cast to ``hmat_dist.tcell``.
    - hmat_dist_cas.tcell: Use CAS to store CT nodes.
    - hmat_dist_casc.tcell: Use CAS to store CT nodes, but in chunks.
4. MPI + OpenMP
    - hmat_div_omp.c: The OpenMP implementation we mentioned in paper.

## 2. matrix partitioning + filling
1. Sequential
    - hmat_filling.c: Sequential version of C which do filling after all leaf-nodes are created.
    - hmat_array_filling_wBCT.c: Sequential version of C which do filling a leaf-node is created.
2. MPI + OpenMP
    - hmat_array_filling_MPI.c: Parallelized ``hmat_filling.c`` using MPI.
    - hmat_array_filling_dynamic.c: Parallelized ``hmat_filling.c`` using MPI and OpenMP with dynamic scheduling.

# Contributors
- [Zhengyang Bai](https://github.com/simon2)
- [Tasuku Hiraishi](https://github.com/tastasgit)
- Hiroshi Nakashima
- Akihiro Ida
- Masahiro Yasugi
- Keiichiro Fukazawa
