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

We also test a hybrid approach that attempts to combine the best properties of the two approaches.
It starts with the linear method where we iterate over vector `X` to find a particular element of vector `Y`.
However, instead of just iterating over consecutive elements in `X`, the steps become exponentially larger.
This is effectively the binary search in reverse, starting from the left-most edge and progressing back to the midpoint.
Once we overstep, we then perform an actual binary search in the relevant interval of `X` to find an exact match to the current element of `Y`.
We then repeat this process for the next element of `Y` until we finish iterating over either of the vectors.
The idea is to achieve linear time when both vectors are highly intermingled,
but to use the exponential step-up and binary search to short-circuit unnecessary iterations when there are gaps in `Y`.

## Results

All timings here are presented in milliseconds on a Dell i7 laptop with 16 GB RAM running Ubuntu 20.04.

In the default parameters, every 10th element (on average) is non-zero, and we want to also extract every 10th element.
This means that the sparse index and extraction vectors are comparable in size, so it is not surprising that the linear method performs better than the binary search.

```console
$ ./build/extractor
Testing a 10000 x 10000 matrix with a density of 0.1
Using a step size of 10 from 0 to 10000
Linear time: 73 for 999712 sum
Binary time: 157 for 999712 sum
Hybrid time: 74 for 999712 sum
```

If we reduce the step size, we increase the size of the extraction vector.
This penalizes the binary search implementation, which is hard-coded to iterate over the (now longer) extraction vector.

```console
$ ./build/extractor --step 1
Testing a 10000 x 10000 matrix with a density of 0.1
Using a step size of 1 from 0 to 10000
Linear time: 185 for 10000475 sum
Binary time: 1034 for 10000475 sum
Hybrid time: 132 for 10000475 sum
```

Conversely, if we increase the step size and the density, we reduce the size of the extraction vector and increase the size of the sparse index vector.
This now favors the binary search.

```console
$ ./build/extractor --density 1 --step 100
Testing a 10000 x 10000 matrix with a density of 1
Using a step size of 100 from 0 to 10000
Linear time: 69 for 1000000 sum
Binary time: 43 for 1000000 sum
Hybrid time: 33 for 1000000 sum
```

Pleasantly enough, the hybrid approach is the best method (or close to it) in all scenarios.
Its superiority over the two pure approaches is somewhat interesting but is probably just a consequence of the implementation;
for linear search, the hybrid approach has hard-coded handling of the "step of 1" scenario,
while for the binary search, the hybrid approach focuses on the left-most interval and avoids redundant traversal of the right-most elements.

Note that there is some opportunity for optimizing all approaches by ensuring that the outer iteration is done on the shorter vector.
This ensures that we process more elements during each invocation of the tight inner loops (or the binary search),
reducing the overhead from the surrounding code.
For example, inverting the `--density 1 --step 100` run above, we get a ~2-fold slow-down when iterating over the larger vector:

```console
$ ./build/extractor --density 0.01 --step 1 
Testing a 10000 x 10000 matrix with a density of 0.01
Using a step size of 1 from 0 to 10000
Linear time: 124 for 998183 sum
Binary time: 644 for 998183 sum
Hybrid time: 77 for 998183 sum
```

Or alternatively, a 2-fold speed-up after inverting the `--step 1` run:

```console
$ ./build/extractor --step 10 --density 1 
Testing a 10000 x 10000 matrix with a density of 1
Using a step size of 10 from 0 to 10000
Linear time: 65 for 10000000 sum
Binary time: 143 for 10000000 sum
Hybrid time: 71 for 10000000 sum
```

## Build instructions

Just use the usual CMake process:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

