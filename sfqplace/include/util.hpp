#ifndef SFQPLACE_UTIL_HPP
#define SFQPLACE_UTIL_HPP

#include <unordered_set>

template<typename T>
std::unordered_set<T> setIntersection(std::unordered_set<T> &set1, std::unordered_set<T> set2) {
    for (const auto &item : set1) {
        if (set2.find(item) != set2.end()) {
            set2.erase(item);
        }
    }

    return set2;
};

#endif // SFQPLACE_UTIL_HPP
