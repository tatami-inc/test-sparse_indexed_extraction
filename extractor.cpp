#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

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

    while (1) {
        auto exval = extract[j];
        auto curval = current[k];

        if (exval < curval) {
            while (1) {
                ++j;
                if (j == num) {
                    return;
                }
                if (extract[j] >= curval) {
                    break;
                }
            }

        } else if (exval > curval) {
            while (1) {
                ++k;
                if (k == end) {
                    return;
                }
                if (exval <= current[k]) {
                    break;
                }
            }

        } else {
            ++collected;
            ++k;
            if (k == end) {
                return;
            }
            ++j;
            if (j == num) {
                return;
            }
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
        auto new_k = std::lower_bound(current.begin() + last_k, current.begin() + k, limit) - current.begin();
        if (new_k < k) {
            if (current[new_k] == limit) {
                ++collected;
                k = new_k + 1;
            } else {
                k = new_k;
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

/** LOOKUP **/

struct LookupTable {
    std::vector<unsigned char> present;
    size_t offset = 0;
};

LookupTable create_lookup_table(const std::vector<int>& extract) {
    LookupTable output;
    if (!extract.empty()) {
        output.offset = extract.front();
        size_t allocation = extract.back() - output.offset + 1;
        output.present.resize(allocation);
        for (auto i : extract) {
            output.present[i - output.offset] = 1;
        }
    }
    return output;
}

size_t collect_lookup(const std::vector<std::vector<int> >& indices, const LookupTable& lookup) {
    size_t collected = 0;
    size_t max = lookup.present.size();
    for (const auto& current : indices) {
        for (auto x : current) {
            // Deliberately creating a branch here, as actual applications will be
            // more complicated than counting the number of discovered elements.
            size_t i = x - lookup.offset;
            if (i < max && lookup.present[i]) {
                ++collected;
            }
        }
    }
    return collected;
}

/** MAIN **/

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

    size_t total_sum = collect_linear(indices, extract);
    std::cout << "Expecting a sum of " << total_sum << std::endl;

    // Running through the possibilities.
    ankerl::nanobench::Bench().run("linear", [&](){
        auto collected = collect_linear(indices, extract);
        if (total_sum != collected) {
            std::cerr << "WARNING: different result from linear access (" << collected << ")" << std::endl;
        }
    });

    ankerl::nanobench::Bench().run("binary", [&](){
        auto collected = collect_pure_binary(indices, extract);
        if (total_sum != collected) {
            std::cerr << "WARNING: different result from binary access (" << collected << ")" << std::endl;
        }
    });

    ankerl::nanobench::Bench().run("hybrid", [&](){
        auto collected = collect_hybrid(indices, extract);
        if (total_sum != collected) {
            std::cerr << "WARNING: different result from hybrid access (" << collected << ")" << std::endl;
        }
    });

    auto tab = create_lookup_table(extract);
    ankerl::nanobench::Bench().run("lookup", [&](){
        auto collected = collect_lookup(indices, tab);
        if (total_sum != collected) {
            std::cerr << "WARNING: different result from lookup access (" << collected << ")" << std::endl;
        }
    });
}
