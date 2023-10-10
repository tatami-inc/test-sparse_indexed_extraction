#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include <cmath>
#include <chrono>
#include <vector>
#include <queue>
#include <random>

size_t collect_linear(const std::vector<std::vector<int> >& indices, const std::vector<int>& extract) {
    size_t num = extract.size();
    if (num == 0) {
        return 0;
    }

    size_t collected = 0;
    for (int c = 0, nc = indices.size(); c < nc; ++c) {
        auto& current = indices[c];

        size_t j = 0;
        if (extract[0]) {
            j = std::lower_bound(current.begin(), current.end(), extract[0]) - current.begin();
        }
        size_t end = current.size();
        size_t k = 0;

        for (; j < end; ++j) {
            auto limit = current[j];
            while (k < num && extract[k] < limit) {
                ++k;
            }
            if (k == num) {
                break;
            }
            if (extract[k] == limit) {
                ++collected;
                ++k;
            }
        }
    }

    return collected;
}

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

size_t intlog2(size_t x) {
    size_t out = 1; // roughly rounding up the log2, ensure that it can't be zero.
    while (x >>= 1) {
        ++out;
    }
    return out;
}

size_t collect_hybrid(const std::vector<std::vector<int> >& indices, const std::vector<int>& extract) {
    size_t num = extract.size();
    if (num == 0) {
        return 0;
    }

    size_t collected = 0;
    for (int c = 0, nc = indices.size(); c < nc; ++c) {
        auto& current = indices[c];
        if (current.empty()) {
            continue;
        }

        // We start by figuring out the intersection of ranges.
        if (current.front() > extract.back() || extract.front() > current.back()) {
            continue;
        }

        auto cur_start = current.begin(), cur_end = current.end();
        auto ext_start = extract.begin(), ext_end = extract.end();
        if (current.front() > extract.front()) {
            ext_start = std::lower_bound(ext_start, ext_end, current.front());
        } else if (current.front() < extract.front()) {
            cur_start = std::lower_bound(cur_start, cur_end, extract.front());
        }

        if (current.back() > extract.back()) {
            ext_end = std::lower_bound(ext_start, ext_end, current.back() + 1);
        } else if (current.back() < extract.back()) {
            cur_end = std::lower_bound(cur_start, cur_end, extract.back() + 1);
        }

        size_t ext_count = ext_end - ext_start;
        size_t cur_count = cur_end - cur_start;

        if (ext_count > cur_count * intlog2(ext_count)) {
            // Doing a binary search for each element in 'current'.
            for (; cur_start != cur_end; ++cur_start) {
                auto curval = *cur_start;
                ext_start = std::lower_bound(ext_start, ext_end, curval);
                if (ext_start == ext_end) {
                    break;
                } else if (*ext_start == curval) {
                    ++collected;
                }
            }

        } else if (cur_count > ext_count * intlog2(cur_count)) {
            // Doing a binary search for each element in 'extract'.
            for (; ext_start != ext_end; ++ext_start) {
                auto extval = *ext_start;
                cur_start = std::lower_bound(cur_start, cur_end, extval);
                if (cur_start == cur_end) {
                    break;
                } else if (*cur_start == extval) {
                    ++collected;
                }
            }

        } else if (ext_count > cur_count) {
            // Inner loop over 'extract', to reduce loop restarts.
            for (; cur_start != cur_end; ++cur_start) {
                auto limit = *cur_start;
                while (ext_start != ext_end && *ext_start < limit) {
                    ++ext_start;
                }
                if (ext_start == ext_end) {
                    break;
                }
                if (*ext_start == limit) {
                    ++collected;
                    ++ext_start;
                }
            }

        } else {
            // Inner loop over 'current', to reduce loop restarts.
            for (; ext_start != ext_end; ++ext_start) {
                auto limit = *ext_start;
                while (cur_start != cur_end && *cur_start < limit) {
                    ++cur_start;
                }
                if (cur_start == cur_end) {
                    break;
                }
                if (*cur_start == limit) {
                    ++collected;
                    ++cur_start;
                }
            }
        }
    }

    return collected;
}


int main(int argc, char* argv []) {
    CLI::App app{"Expanded testing checks"};
    double density;
    app.add_option("-d,--density", density, "Density of the expanded sparse matrix")->default_val(0.1);
    int nr;
    app.add_option("-r,--nrow", nr, "Number of rows")->default_val(10000);
    int nc;
    app.add_option("-c,--ncol", nc, "Number of columns")->default_val(10000);
    double start;
    app.add_option("--start", start, "Start of the extraction, as a fraction of the number of rows")->default_val(0.1);
    double end;
    app.add_option("--end", end, "End of the extraction, as a fraction of the number of rows")->default_val(0.9);
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

    // Running through the possibilities.
    {
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_linear(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << "Boosted linear time: " << duration.count() << " for " << collected << " sum" << std::endl;
    }

    {
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_pure_binary(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << "Pure binary time: " << duration.count() << " for " << collected << " sum" << std::endl;
    }

    {
        auto tstart = std::chrono::high_resolution_clock::now();
        auto collected = collect_hybrid(indices, extract);
        auto tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << "Hybrid time: " << duration.count() << " for " << collected << " sum" << std::endl;
    }
}
