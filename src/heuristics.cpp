#include "heuristics.hpp"
#include "librarySPP.hpp"

#include <cmath>
#include <omp.h>

std::tuple<char*, int, char*> GreedyRandomized(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U,
        const float alpha) {
    bool valid;
    int i(0), j(0), k(0), s(0), e(0), min_u, max_u;
    float limit(0.0f);
    std::vector<int> RCL;
    char *x = new char[n], *column = new char[m];
    for(i = 0; i < n; i++) x[i] = 0;
    for(j = 0; j < m; j++) column[j] = 0;

    int* u_order = argsort(n, U); // DON'T FORGET TO DELETE

    k = 0;
    while(s != m && k < n) {
        // Indices of max and min utilities
        for(j = 0, max_u = 0; j < n && u_order[max_u] == -1; j++)
            max_u = j;
        for(j = n-1, min_u = 0; j >= 0 && u_order[min_u] == -1; j--)
            min_u = j;
        limit = U[u_order[min_u]] + alpha * (U[u_order[max_u]] - U[u_order[min_u]]);
        for(j = 0; j < n; j++)
            // If variable's index is candidate and the utility
            // greater than limit then we add the index in RCL
            if(u_order[j] != -1 && U[u_order[j]] >= limit)
                RCL.push_back(j);

        // Select an element e from RCL at random
        e = (RCL.size()) ? *select_randomly(RCL.begin(), RCL.end()) : max_u;
        for(j = 0, valid = true; j < m && valid; j++)
            valid = column[j] + A[INDEX(u_order[e], j)] <= 1;
        for(j = 0, s = 0; valid && j < m; s += column[j], j++)
            column[j] += A[INDEX(u_order[e], j)];
        x[u_order[e]] = valid, u_order[e] = -1;
        k += 1; RCL.clear();
    }

    delete[] u_order;
    return std::make_tuple(x, dot(n, x, C), column);
}

void GreedyImprovement(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column) {
    int i(2);
    bool (*f[3])(int, int, const int*, const char*, char*, int*, bool, char*) = {
            zero_oneExchange,
            one_oneExchange,
            two_oneExchange
        };

    // We modify x and z directly (no copy)
    while(i >= 0){
        if(!f[i](m, n, C, A, x, z, deep, column)) i--;
    }
}

void ReactiveGRASP(
        const int m,
        const int n,
        const int* C,
        const char* A,
        const float* U,
        std::vector<int>& zInits,
        std::vector<int>& zAmels,
        std::vector<int>& zBests,
        const std::vector<double>& alpha,
        std::vector<double>& proba,
        const int probaUpdate,
        const double delta,
        const int nbIter,
        const bool deep,
        const bool parallel) {
    int iter(0), zBest(-1), update(0);
    std::vector<std::vector<int>> pool(alpha.size());
    std::vector<double> valuation(alpha.size(), 0);

    #pragma omp parallel for ordered if(parallel)
    for(iter = 0; iter < nbIter; iter++) {
        char *x(nullptr), *column(nullptr);
        int i(0), j(0);
        double sel_alpha(-1), idx((double)rand() / RAND_MAX), s(0);
        for(i = 0; i < (int)proba.size() && sel_alpha == -1; i++) {
            s += proba[i];
            if(idx < s) sel_alpha = proba[i];
        }
        if(i == (int)proba.size()) i--;
        std::tie(x, zInits[iter], column) = GreedyRandomized(m, n, C, A, U, sel_alpha);
        zAmels[iter] = zInits[iter];
        GreedyImprovement(m, n, C, A,
                x, &zAmels[iter], deep, column);

        // Section de code difficilement parallélisable
        if(update == probaUpdate) {
            #pragma omp critical
            {
                int zMax = *std::max_element(zAmels.begin(), zAmels.end()),
                    zMin = *std::min_element(zAmels.begin(), zAmels.end());
                for(j = 0; j < (int)alpha.size(); j++) {
                    double mean(0); for(auto e : pool[j]) mean += e;
                    mean = (pool[j].size() ? mean/pool[j].size() : 0);
                    valuation[j] = std::pow((mean - zMin)/(zMax-zMin), delta);
                    pool[j].clear();
                }
                double sum(0); for(auto e : valuation) sum += e;
                for(j = 0; j < (int)proba.size(); j++) {
                    proba[j] = valuation[j]/sum;
                    proba[j] = (std::isnan(proba[j]) ? 0.0 : proba[j]);
                }
                update = 0;
            }
        }

        pool[i].push_back(zAmels[iter]), update++;

        /* MOST IMPORTANT SECTION */
        if(x) delete[] x, x = nullptr;
        if(column) delete[] column, column = nullptr;
    }

    // Compute zBests using zAmels
    for(iter = 0; iter < nbIter; iter++) {
        zBest = std::max(zBest, zAmels[iter]);
        zBests[iter] = zBest;
    }
}
