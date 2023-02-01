# Synopsis

H-div is a parallelized implementation of matrix partitioning in construction of hierarchical matrices using Cilk Plus and Tascell based on a sequential Fortran implementation in HACApK library coded by Prof. Akihito Ida.

To learn more about HACApK, please visit [HACApK introduction](http://ppopenhpc.cc.u-tokyo.ac.jp/ppopenhpc/2017/01/31/ppopen-applbem-ver-0-5-0/) or [HACApK github repository](https://github.com/Post-Peta-Crest/ppOpenHPC/tree/MATH/HACApK).

To learn more about Tascell, please visit [Tascell introduction](http://ais.sys.i.kyoto-u.ac.jp/~task/tascell/index.html) or [Tascell github repository](https://github.com/tascell/sc-tascell).

For details of this implementation, please read paper [_Parallelization of Matrix Partitioning in Construction of Hierarchical Matrices using Task Parallel Languages_](https://www.jstage.jst.go.jp/article/ipsjjip/27/0/27_840/_article/-char/ja/)

# Requirements
- An Intel multi-core CPU
- To run Cilk Plus versions, Intel C++ Compiler version >= 17
- To run tascell versions, 
  - tascell compiler version later than Jan 21, 2019.
  - GCC version >= 4.8.5

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
- hmat_div.c: The first sequential implementation using C directly translated from Fortran implementation in HACApK.
- hmat_div.cpp: Sequential implementation of C++.
- hmat_div.sc: Sequential implementation using S-expression-based syntax. This is the base of Tascell.
- hmat_div.tcell: The final implementation using Tascell.
- hmat_div_BCT_cilk_list_reducer.cpp: Tried CILK_LIST_REDUCER, but not included in paper due to bad performance.
- hmat_div_BCT_cilk_malloc.c: Create private BCT array for each worker by using malloc function. Not included in paper due to bad performance.
- hmat_div_CT_cilk_parByLevel.c: Tried to swith parallel code to sequential code by tree level information. Not included in paper due to bad performance.
- hmat_div_array.c: CT is not in linked-tree manner, but use a pre-allocated array. (will be discribed in next paper)
- hmat_div_cilk.c: Final implementation using Cilk Plus.
- hmat_div_direct.c: Sequential C implementation but exchange data elements directly instead of exchange their index. This is the baseline of sequential implementation in paper.
- hmat_div_locality.tcell: Tried to make upper tree levels execute sequentially and execute in parallel in lower levels, for better data locality. However, the result of speedup is not essential.
- hmat_div_omp.c: The OpenMP implementation we mentioned in paper.
# Contributors
- [Zhengyang Bai](https://github.com/simon2)
- [Tasuku Hiraishi](https://github.com/tastasgit)
- Hiroshi Nakashima
- Akihiro Ida
- Masahiro Yasugi
- Keiichiro Fukazawa
