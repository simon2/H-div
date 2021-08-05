# Synopsis

H-div is a parallelized implementation of matrix partitioning in construction of hierarchical matrices using Cilk Plus and Tascell based on a sequential Fortran implementation in HACApK library coded by Prof. Akihito Ida.

To learn more about HACApK, please visit [HACApK introduction](http://ppopenhpc.cc.u-tokyo.ac.jp/ppopenhpc/2017/01/31/ppopen-applbem-ver-0-5-0/).

To learn more about Tascell, please visit [Tascell introduction](http://ais.sys.i.kyoto-u.ac.jp/~task/tascell/index.html) or [Tascell github repository](https://github.com/tascell/sc-tascell).

For details of this implementation, please read paper [_Parallelization of Matrix Partitioning in Construction of Hierarchical Matrices using Task Parallel Languages_](https://www.jstage.jst.go.jp/article/ipsjjip/27/0/27_840/_article/-char/ja/)

If you want to cite this paper, please use:
```
@article{Zhengyang Bai2019,
  title={Parallelization of Matrix Partitioning in Construction of Hierarchical Matrices using Task Parallel Languages},
  author={Zhengyang Bai and Tasuku Hiraishi and Hiroshi Nakashima and Akihiro Ida and Masahiro Yasugi},
  journal={Journal of Information Processing},
  volume={27},
  number={ },
  pages={840-851},
  year={2019},
  doi={10.2197/ipsjjip.27.840}
}
```
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
# Contributors
- [Zhengyang Bai](https://github.com/simon2)
- [Tasuku Hiraishi](https://github.com/tastasgit)
- Hiroshi Nakashima
- Akihiro Ida
- Masahiro Yasugi
