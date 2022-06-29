#pragma once
#include <iostream>
#include <vector>

class MCEngine {
public:
    template <class Function, class T>
    std::vector<T> Simulate(Function func, uint64_t n_simulations = 10000u) {
        std::vector<T> data;

        for (auto i = 0u; i < n_simulations; ++i) {
            data.emplace_back(func());
        }

        return data;
    }
};
