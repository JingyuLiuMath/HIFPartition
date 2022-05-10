# HIFPartition

HIFPartition is a fast solver for SPD (symmetric positive definite) linear systems. The MF algorithm computes the LDL decomposition of a SPD matrix exactly. The HIF algorithm computes the approximate generalized LDL decomposition of a matrix but is faster than MF. In this code, we do HIF via a graph partition given by [metismex](https://github.com/YingzhouLi/metismex) or [meshpart](https://github.com/YingzhouLi/meshpart).

## Authors

* Jingyu Liu, Fudan University, 381258337@qq.com

## Aknowledgement

Many thanks to [Yingzhou Li](https://www.yingzhouli.com/). Without his guidance and help, the code can't be finished.
