#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include <cmath>
#include <chrono>
#include <vector>
#include <queue>
#include <random>

/** LINEAR **/

void collect_linear_internal(const std::vector<int>& current, const std::vector<int>& extract, size_t& collected) {
    size_t num = extract.size();
    size_t j = 0;
    if (current[0]) {
        j = std::lower_bound(extract.begin(), extract.end(), current[0]) - extract.begin();
    }

    size_t end = current.size();
    size_t k = 0;
    if (extract[0]) {
        k = std::lower_bound(current.begin(), current.end(), extract[0]) - current.begin();
    }

    for (; j < num; ++j) {
        auto limit = extract[j];

        // Skip the loop if the current element is already a match or exceeds the limit.
        if (limit < current[k]) {
            continue;
        } else if (limit == current[k]) {
            ++k;
            ++collected;
            if (k == end) {
                return;
            }
            continue;
        }

        while (1) {
            ++k;
            if (k == end) {
                return; // an *_internal function allows u to break out of the outer loop immediately by just returning.
            } 
            if (current[k] >= limit) {
                break;
            }
        }

        if (current[k] == limit) {
            ++collected;
            ++k;
        }
    }
}

size_t collect_linear(const std::vector<std::vector<int> >& indices, const std::vector<int>& extract) {
    size_t num = extract.size();
    if (num == 0) {
        return 0;
    }

    size_t collected = 0;
    for (int c = 0, nc = indices.size(); c < nc; ++c) {
        collect_linear_internal(indices[c], extract, collected);
    }

    return collected;
}

/** BINARY **/

size_t collect_pure_binary(const std::vector<std::vector<int> >& indices, const std::vector<int>& extract) {
    if (extract.empty()) {
        return 0;
    }

    size_t collected = 0;
    for (int c = 0, nc = indices.size(); c < nc; ++c) {
        auto& current = indices[c];
        auto sofar = current.begin(), end = current.end();
        for (auto x : extract) {
            sofar = std::lower_bound(sofar, end, x);
            if (sofar == end) {
                break;
            } else if (*sofar == x) {
                ++collected;
            }
        }
    }

    return collected;
}

/** HYBRID **/

void collect_hybrid_internal(const std::vector<int>& current, const std::vector<int>& extract, size_t& collected) {
    size_t num = extract.size();
    size_t j = 0;
    if (current[0]) {
        j = std::lower_bound(extract.begin(), extract.end(), current[0]) - extract.begin();
    }

    size_t end = current.size();
    size_t k = 0;
    if (extract[0]) {
        k = std::lower_bound(current.begin(), current.end(), extract[0]) - current.begin();
    }

    for (; j < num; ++j) {
        auto limit = extract[j];

        // Handle the common case of the current[k] already exceeding/equalling the limit.
        if (current[k] > limit) {
            continue;
        } else if (current[k] == limit) {
            ++collected;
            ++k;
            if (k == end) {
                return;
            }
            continue;
        }

        // Use an exponential step-up, starting with +1, then +2, then +4,
        // and so on.  This could be interpreted as the reverse of a binary
        // search that terminates at the left-most edge. We special-case
        // the initial step of +1 as it's pretty common.
        ++k;
        if (k == end) {
            return;
        }
        if (current[k] > limit) {
            continue;
        } else if (current[k] == limit) {
            ++collected;
            ++k;
            if (k == end) {
                return;
            }
            continue;
        } 

        size_t step = 1, last_k = k;
        do {
            step <<= 1; // i.e., step of 2, then 4, then 8 ... 
            if (step >= end - k) { // avoid issues with overflow.
                k = end;
                break;
            }
            last_k = k;
            k += step;
        } while (current[k] < limit);

        if (k < end && current[k] == limit) {
            ++collected;
            ++k;
            if (k == end) {
                return;
            } 
            continue;
        }

        // Perform a binary search to trim down any overshooting after the
        // step-up. If a binary search is treated as a decision tree, we
        // basically just walked up the tree from the left-most edge (i.e.,
        // the 'k' at the start) to some intermediate node (or the root)
        // and now we're walking back down to find the 'limit'.
        size_t right = k;
        k = last_k + 1; // 'current[last_k]' must be less than 'limit', as we would have otherwise escaped the loop sooner.

        while (k < right) {
            size_t mid = k + ((right - k) >> 1); 
            auto midval = current[mid];
            if (midval == limit) {
                ++collected;
                k = mid + 1;
                break;
            } else if (midval > limit) {
                right = mid;
            } else {
                k = mid + 1;
            }
        }

        if (k == end) {
            return;
        }
    }
}

size_t collect_hybrid(const std::vector<std::vector<int> >& indices, const std::vector<int>& extract) {
    size_t num = extract.size();
    if (num == 0) {
        return 0;
    }

    size_t collected = 0;
    for (int c = 0, nc = indices.size(); c < nc; ++c) {
        collect_hybrid_internal(indices[c], extract, collected);
    }

    return collected;
}


int main(int argc, char* argv []) {
    CLI::App app{"Expanded testing checks"};
    double density;
    app.add_option("-d,--density", density, "Density of the expanded sparse matrix")->default_val(0.1);
    int nr;
    app.add_option("-r,--nrow", nr, "Number of rows")->default_val(50000);
    int nc;
    app.add_option("-c,--ncol", nc, "Number of columns")->default_val(10000);
    double start;
    app.add_option("--start", start, "Start of the extraction, as a fraction of the number of rows")->default_val(0);
    double end;
    app.add_option("--end", end, "End of the extraction, as a fraction of the number of rows")->default_val(1);
    int step;
    app.add_option("--step", step, "Step size of the extraction, in terms of number of rows")->default_val(10);
    CLI11_PARSE(app, argc, argv);

    std::cout << "Testing a " << nr << " x " << nc << " matrix with a density of " << density << std::endl;

    // Simulating a set of sparse vectors.
    std::mt19937_64 generator(1234567);
    std::uniform_real_distribution<double> distu;

    std::vector<std::vector<int> > indices(nc);
    for (int c = 0; c < nc; ++c) {
        auto& current = indices[c];
        for (int r = 0; r < nr; ++r) {
            if (distu(generator) <= density) {
                current.push_back(r);
            }
        }
    }

    // Simulating the queries.
    std::vector<int> extract;
    int true_start = start * nr;
    int true_end = end * nr;
    for (int r = true_start; r < true_end; r += step) {
        extract.push_back(r);
    }
    std::cout << "Using a step size of " << step << " from " << true_start << " to " << true_end << std::endl;

    // Running through the possibilities.
    {
        std::cout << "Linear time: ";
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_linear(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << duration.count() << " for " << collected << " sum" << std::endl;
    }

    {
        std::cout << "Binary time: ";
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_pure_binary(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << duration.count() << " for " << collected << " sum" << std::endl;
    }

    {
        std::cout << "Hybrid time: ";
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_hybrid(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << duration.count() << " for " << collected << " sum" << std::endl;
    }
}
