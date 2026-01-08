//
// Created by Alejandro Zeise on 12/12/23.
//

#ifndef PA3ANALYTICPLACEMENT_UTIL_HPP
#define PA3ANALYTICPLACEMENT_UTIL_HPP

#include <utility>
#include <cstddef>
#include <functional>

struct firstElementOfPairComparator
{
    bool operator()(const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) const
    {
        return lhs.first < rhs.first;
    }
};

// https://stackoverflow.com/questions/20590656/error-for-hash-function-of-pair-of-ints
struct pairHashInteger final {
    size_t operator()(const std::pair<int, int>& p) const noexcept {
        size_t hash = std::hash<int>{}(p.first);
        hash <<= sizeof(size_t) * 4;
        hash ^= std::hash<int>{}(p.second);
        return std::hash<size_t>{}(hash);
    }
};

#endif //PA3ANALYTICPLACEMENT_UTIL_HPP
