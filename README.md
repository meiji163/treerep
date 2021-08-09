# Treerep
C++ implementation of Sonathalia & Gilbert "*Tree! I am no Tree! I am a Low Dimensional Hyperbolic Embedding,"* <sup>[1](#1)</sup> with a thin python 3 wrapper. 

Treerep takes an N-point metric (matrix of pairwise distances) and computes a weighted tree that approximates it. You can use this tree directly, or embed it in hyperbolic space<sup>[2](#2)</sup>.
If the original metric is [hyperbolic](https://en.wikipedia.org/wiki/Hyperbolic_metric_space) (roughly "tree-like") then Treerep produces an embedding with low distortion. Recent research<sup>[3](#3)</sup> <sup>[4](#4)</sup> <sup>[5](#5)</sup> shows hyperbolic embeddings outperform Euclidean embeddings ones for many types of hierarchical data.

## Install
The C++11 source can be used with no external dependencies. The python API requires pybind11, which can be downloaded by cloning the submodule
```shell
$ git clone --recursive https://github.com/meiji163/treerep.git
```
then compile with cmake
```shell
$ cd treerep && mkdir build && cd build
$ cmake ..
$ make
$ make install #if you want to install to your python lib
```

## Usage
To save space, the N-point metric is represented as a length N*(N-1)/2 array (see [scipy pdist](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist)). Input can be a python list or numpy array.

Here is Sarich's "immunological distance" example from Sonthalia & Gilbert's paper.

```python
from trep import *
from trep_utils import load_mat
metric, N = load_mat("../data/immune.mtx")
print(metric)
# array([ 32.,  48.,  51.,  50.,  48.,  98., 148.,  26.,  34.,  29.,  33.,
#        84., 136.,  42.,  44.,  44.,  92., 152.,  44.,  38.,  86., 142.,
#        42.,  89., 142.,  90., 142., 148.])

tree, weights = treerep(metric,N)
print(weights) 
# {(0, 10): 24.0, (1, 10): 8.0, (2, 11): 19.25, (3, 9): 18.25, (4, 12): 21.25, 
#  (5, 8): 17.5, (6, 9): 67.75, (7, 13): 121.9375, (8, 9): 2.25, (8, 13): 1.0625, 
#  (10, 11): 1.75, (11, 12): 2.0, (12, 13): 1.6875}
```
the tree looks like this (edges not to scale)

![sarich tree](https://i.imgur.com/XsCesv1.png)

Treerep is a randomized algorithm, so you can produce multiple trees and choose the best.

## Todo
* multithread sorting
* distortion/ MAP

## References
<div><a name="1">1</a>: https://arxiv.org/abs/2005.03847
<div><a name="2">2</a>: https://arxiv.org/abs/1804.03329
<div><a name="3">3</a>: https://arxiv.org/abs/1904.02239
<div><a name="4">4</a>: https://arxiv.org/abs/1705.08039
<div><a name="5">5</a>: https://arxiv.org/abs/1810.06546

