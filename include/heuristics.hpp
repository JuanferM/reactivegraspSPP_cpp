#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "movements.hpp"

#include <string>
#include <cstring>

// Greedy randomized construction of a feasible solution
std::tuple<char*, int, char*> GreedyRandomized(
    int m,
    int n,
    const int* C,
    const char* A,
    const float* U,
    const float alpha);

// Greedy improvement of a feasible solution through (deep) local search
void GreedyImprovement(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int z,
    bool deep = true,
    char* column = nullptr);

// GRASP for the Set Packing Problem
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
    int nbIter = 100,
    bool deep = true,
    bool parallel = true);

#endif /* end of include guard: HEURISTICS_H */
