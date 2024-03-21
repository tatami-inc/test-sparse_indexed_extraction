# Performance testing for indexed sparse extraction

## Overview

When dealing with sparse matrices in **tatami**,
a common task is to extract elements from one sorted vector (the sparse indices) based on an another sorted vector (the extraction indices).
There are several broad strategies to do so:

- **Linear**: we scan through both vectors simultaneously, comparing elements at each step.
  This takes advantage of the fact that both vectors are sorted and runs in linear time with respect to the larger vector.
- **Binary**: we iterate over one of the vectors, performing a binary search on the other vector.
  This runs in linear time with respect to the iterated vector and is `O(log(n))` with respect to the other vector.
- **Hybrid**: this attempts to combine the best properties of the linear and binary approaches.
  It starts with the linear method where we iterate over vector `X` to find a particular element of vector `Y`.
  However, instead of just iterating over consecutive elements in `X`, the steps become exponentially larger.
  This is effectively the binary search in reverse, starting from the left-most edge and progressing back to the midpoint.
  Once we overstep, we then perform an actual binary search in the relevant interval of `X` to find an exact match to the current element of `Y`.
  We then repeat this process for the next element of `Y` until we finish iterating over either of the vectors.
  The idea is to achieve linear time when both vectors are highly intermingled,
  but to use the exponential step-up and binary search to short-circuit unnecessary iterations when there are gaps in `Y`.
- **Lookup**: we use one vector to build a lookup table of booleans indicating whether a value is present.
  We then scan through the other vector, inspecting the lookup table for each element to identify a match.
  This is linear with respect to the iterated vector, scaled by the cost of each lookup.
  The lookup table construction is linear with respect to its input vector,
  but as the table is created once and re-used for multiple matrix rows/columns, 
  its construction time can essentially be ignored.
  We'll use a simple `std::vector` as the lookup table, favoring speed at the cost of some memory efficiency.

This repository builds an executable that compares these approaches under various scenarios.

## Results

All timings here are presented in milliseconds on a Dell i7 laptop with 16 GB RAM running Ubuntu 22.04, 
compiled with GCC 11.4.0 in release mode (`-O3`).

With the default parameters, every 10th element (on average) is non-zero, and we want to also extract every 10th element.
The lookup approach is the fastest by an order of magnitude, emphasizing the benefits of spending some memory on a lookup table.
The linear method performs better than the binary search as the sparse index and extraction vectors are comparable in size.

```console
$ ./build/extractor
Testing a 50000 x 10000 matrix with a density of 0.1
Using a step size of 10 from 0 to 50000
Expecting a sum of 4996708

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      367,654,826.00 |                2.72 |    0.7% |      4.05 | `linear`
|      655,226,331.00 |                1.53 |    1.3% |      7.20 | `binary`
|      363,575,712.00 |                2.75 |    0.2% |      4.00 | `hybrid`
|       30,349,434.00 |               32.95 |    0.1% |      0.34 | `lookup`
```

If we reduce the step size, we increase the size of the extraction vector.
This penalizes the binary search implementation, which is hard-coded to iterate over the (now longer) extraction vector.

```console
$ ./build/extractor --step 1
Testing a 50000 x 10000 matrix with a density of 0.1
Using a step size of 1 from 0 to 50000
Expecting a sum of 49995342

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      688,490,196.00 |                1.45 |    1.2% |      7.55 | `linear`
|    3,850,187,517.00 |                0.26 |    0.5% |     42.33 | `binary`
|      611,670,891.00 |                1.63 |    0.6% |      6.73 | `hybrid`
|       31,169,766.00 |               32.08 |    0.7% |      0.35 | `lookup`
```

Conversely, if we increase the step size and the density, we reduce the size of the extraction vector and increase the size of the sparse index vector.
Now binary search is starting to be less bad:

```console
$ ./build/extractor --density 0.3 --step 20
Testing a 50000 x 10000 matrix with a density of 0.3
Using a step size of 20 from 0 to 50000
Expecting a sum of 7493525

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      339,974,770.00 |                2.94 |    0.3% |      3.76 | `linear`
|      509,600,720.00 |                1.96 |    0.5% |      5.63 | `binary`
|      402,051,880.00 |                2.49 |    0.1% |      4.42 | `hybrid`
|       96,752,524.00 |               10.34 |    2.9% |      1.02 | `lookup`
```

Or if we take it to the extreme, we can favor the binary search:

```console
$ ./build/extractor --density 1 --step 100
Testing a 50000 x 10000 matrix with a density of 1
Using a step size of 100 from 0 to 50000
Expecting a sum of 5000000

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      421,827,513.00 |                2.37 |    0.5% |      4.65 | `linear`
|      224,899,200.00 |                4.45 |    0.9% |      2.46 | `binary`
|      130,942,224.00 |                7.64 |    1.1% |      1.44 | `hybrid`
|      291,717,233.00 |                3.43 |    0.5% |      3.22 | `lookup`
```

Pleasantly enough, the hybrid approach performs close to (or is) the best method in all scenarios.
Its superiority over the binary search in the last scenario is bcause the hybrid approach focuses on the left-most interval,
thus avoiding redundant traversal of the right-most elements.

## Concluding remarks

If one is willing to spend some memory, the lookup table wins convincingly in most scenarios.
This strategy is suboptimal for dense matrices, but storing dense data in sparse format is already inefficient,
so some minor performance degradation for indexed extraction is the least of our concerns here.
Importantly, the lookup table is very easy to implement, which makes it appealing for general usage in **tatami**.

If we're not willing to spend memory (e.g., if the dimension extent is very large), 
then the linear search performs pretty well, even in the worst-case scenarios that should favor the binary search.
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
