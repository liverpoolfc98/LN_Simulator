#pragma once
#include <cmath>
#include <random>

class MetropolisHasting {
private:
    double m_temperature;
    uint64_t m_internal_iterations_num;
    uint64_t m_external_iterations_num;
public:
    MetropolisHasting(double a_temperature, uint64_t a_internal_iterations_num = 100u,
        uint64_t a_external_iterations_num = 1000u)
        : m_temperature(a_temperature),
          m_internal_iterations_num(a_internal_iterations_num),
          m_external_iterations_num(a_external_iterations_num) {}

    template <class Function, class Generator, class T>
    T MCMC(Function func, Generator gen, T point) const {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> dist(0, 1);

        for (auto i = 0u; i < m_internal_iterations_num; ++i) {
            T transition = gen(point);
            double difference = func(transition) - func(point);

            if (dist(rng) < std::exp(-1.0 / m_temperature * difference)) {
                point = transition;
            }
        }

        return point;
    }

    void ChangeTemperature() {
        m_temperature *= 0.95;
    }

    template <class Function, class Generator, class T, bool log = false>
    T Solver(Function func, Generator gen) {

        T point = gen(T());

        for (auto i = 0u; i < m_external_iterations_num; ++i) {
            auto transition = MCMC(func, gen, point);
            if (transition == point) {
                return point;
            }
            point = transition;

            if (log) {
                std::cout << "Point: " << point << std::endl << "Stage: " << i << std::endl
                    << "Functional: " << func(point) << std::endl << std::endl;
            }

            ChangeTemperature();
        }

        return point;
    }
};
