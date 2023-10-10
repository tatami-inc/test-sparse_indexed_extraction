# Performance testing for indexed sparse extraction

## Overview

When dealing with sparse matrices in **tatami**,
a common task is to extract elements from one sorted vector (the sparse indices) based on an another sorted vector (the extraction indices).
There are two broad strategies to do so:

- **Linear**: we scan through both vectors simultaneously, comparing elements at each step.
  This takes advantage of the fact that both vectors are sorted and runs in linear time with respect to the larger vector.
- **Binary**: we iterate over one of the vectors, performing a binary search on the other vector.
  This runs in linear time with respect to the iterated vector and is `O(log(n))` with respect to the other vector.

This repository builds an executable that compares these two approaches under various scenarios.
We also test a hybrid approach that switches between the linear and binary methods based on the relative lengths of the two vectors.

## Results

All timings here are presented in milliseconds on a Dell i7 laptop with 16 GB RAM running Ubuntu 20.04.

In the default parameters, every 10th element (on average) is non-zero, and we want to also extract every 10th element.
This means that the sparse index and extraction vectors are comparable in size, so it is not surprising that the linear method performs best.

```console
$ ./build/extractor
Testing a 10000 x 10000 matrix with a density of 0.1
Using a step size of 10 from 0 to 10000
Linear time: 79 for 999712 sum
Binary time: 128 for 999712 sum
Hybrid time: 76 for 999712 sum
```

If we reduce the step size, we increase the size of the extraction vector.
This penalizes the binary search implementation, which is hard-coded to iterate over the (now longer) extraction vector.

```console
$ ./build/extractor --step 1
Testing a 10000 x 10000 matrix with a density of 0.1
Using a step size of 1 from 0 to 10000
Linear time: 134 for 10000475 sum
Binary time: 659 for 10000475 sum
Hybrid time: 130 for 10000475 sum
```

Conversely, if we increase the step size and the density, we reduce the size of the extraction vector and increase the size of the sparse index vector.
This now favors the binary search. 

```console
$ ./build/extractor --density 0.8 --step 50
Testing a 10000 x 10000 matrix with a density of 0.8
Using a step size of 50 from 0 to 10000
Linear time: 121 for 1599921 sum
Binary time: 63 for 1599921 sum
Hybrid time: 68 for 1599921 sum
```

Pleasantly enough, the hybrid approach is or comparable to the best method in all scenarios.

## Comments on the hybrid approach

The logic here is quite simple.
If one vector is considerably longer than the other, we iterate over the latter and perform a binary search on the former.
This allows us to mitigate the longer vector by reducing it to `O(log(n))` time complexity.
On the other hand, if the vectors are mostly of the same length, we perform a linear search
as we are not being strongly penalized by a difference in vector length.

Here, the key question is what is "considerably longer".
Assuming a logarithmic time complexity on average, we switch to the binary search when `m > n * log2(m)` where `m` is the longer vector and `n` is the shorter vector.
The approximate `log2` can be computed relatively quickly for integers with some bit shifting.
In practice, we're probably missing an extra scaling factor here, but that may be implementation- and machine-dependent, so whatever.

We also need to do some work to actually account for the relevant length of each vector.
If the vectors don't actually overlap in terms of their indices, then the relevant length of each vector is zero.
So we always need to execute up to two binary searches to find the boundaries of one vector within the other vector;
we then consider the length of each vector with respect to the number of elements inside the shared interval.

## Build instructions

Just use the usual CMake process:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

