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

All timings here are presented in milliseconds on a Dell i7 laptop with 16 GB RAM running Ubuntu 20.04, compiled with GCC 9.4.0.

In the default parameters, every 10th element (on average) is non-zero, and we want to also extract every 10th element.
This means that the sparse index and extraction vectors are comparable in size, so it is not surprising that the linear method performs better than the binary search.

```console
$ ./build/extractor
Testing a 50000 x 10000 matrix with a density of 0.1
vUsing a step size of 10 from 0 to 50000
Linear time: 321 for 4996708 sum
Binary time: 651 for 4996708 sum
Hybrid time: 354 for 4996708 sum
```

If we reduce the step size, we increase the size of the extraction vector.
This penalizes the binary search implementation, which is hard-coded to iterate over the (now longer) extraction vector.

```console
$ ./build/extractor --step 1
Testing a 50000 x 10000 matrix with a density of 0.1
Using a step size of 1 from 0 to 50000
Linear time: 579 for 49995342 sum
Binary time: 6337 for 49995342 sum
Hybrid time: 601 for 49995342 sum
```

Conversely, if we increase the step size and the density, we reduce the size of the extraction vector and increase the size of the sparse index vector.
Now binary search is starting to be less bad:

```console
$ ./build/extractor --density 0.3 --step 20
Testing a 50000 x 10000 matrix with a density of 0.3
Using a step size of 20 from 0 to 50000
Linear time: 293 for 7493525 sum
Binary time: 595 for 7493525 sum
Hybrid time: 384 for 7493525 sum
```

Or if we take it to the extreme, we can favor the binary search:

```console
$ ./build/extractor --density 1 --step 100
Testing a 50000 x 10000 matrix with a density of 1
Using a step size of 100 from 0 to 50000
Linear time: 253 for 5000000 sum
Binary time: 225 for 5000000 sum
Hybrid time: 132 for 5000000 sum
```

Pleasantly enough, the hybrid approach performs close to (or is) the best method in all scenarios.
Its superiority over the binary search in the last scenario is bcause the hybrid approach focuses on the left-most interval,
thus avoiding redundant traversal of the right-most elements.

Note that we can easily improve the binary search by ensuring that the outer iteration is done on the shorter vector,
so that the logarithmic time complexity can be applied to the larger vector.
For example, we see a major speed-up after inverting the `--step 1 (--density 0.1)` run, though it is still slower than the linear search.

```console
$ ./build/extractor --step 10 --density 1
Testing a 50000 x 10000 matrix with a density of 1
Using a step size of 10 from 0 to 50000
Linear time: 236 for 50000000 sum
Binary time: 790 for 50000000 sum
Hybrid time: 357 for 50000000 sum
```

## Concluding remarks

Despite its simpliciy, the linear search performs pretty well, even in the worst-case scenarios that should favor the binary search.
I assume that this is because branch prediction allows uninteresting elements to be skipped efficiently in the linear method,
whereas the binary search suffers from unpredictable branches.
I'd also speculate that both methods are subject to the same memory bandwidth limits - 
the linear search obviously needs to look at each element,
while the binary search probably operates within a single cache line (and thus does not really skip loading of any element).

The hybrid approach is also good and I suppose we could just use it all the time.
It performs better than the linear approach in a few scenarios and is at least competitive in the other cases.
However, it is a lot more complicated to implement.

## Build instructions

Just use the usual CMake process:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

