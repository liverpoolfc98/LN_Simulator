#pragma once
#include <cmath>

class Market {
private:
    double m_interest_rate;
    double m_current_time;
public:
    Market(double a_interest_rate, double a_current_time = 0.0)
    : m_interest_rate(a_interest_rate), m_current_time(a_current_time) {}

    void Update(double time_interval) {
        m_current_time += time_interval;
    }

    double GetInterestRate() const {
        return m_interest_rate;
    }

    double GetCurrentTime() const {
        return m_current_time;
    }

    double Discount(double value, double initial_time = 0.0) const {
        double discount_factor = std::exp(-(m_current_time - initial_time) * m_interest_rate);

        return value * discount_factor;
    }

    double BondValue(double value, double time_interval) const {
        double inverse_discount_factor = std::exp(time_interval * m_interest_rate);

        return value * inverse_discount_factor;
    }
};