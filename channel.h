#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include <string>

#include "market.h"

class Channel {
private:
    int64_t m_balance_lhs, m_balance_rhs;
    double m_open_time;
    Market* m_market;
public:
    Channel(int64_t a_balance_lhs, int64_t a_balance_rhs, Market* a_market)
    : m_balance_lhs(a_balance_lhs),
      m_balance_rhs(a_balance_rhs),
      m_open_time(0.0),
      m_market(a_market) {}

    bool ProcessPayment(int64_t value) {
        if ((value > 0 && m_balance_lhs < value) || (value < 0 && m_balance_rhs < -value)) {
            return false;
        }
        m_balance_lhs -= value, m_balance_rhs += value;
        return true;
    }

    void Reset(int64_t a_balance_lhs, int64_t a_balance_rhs) {
        m_balance_lhs = a_balance_lhs, m_balance_rhs = a_balance_rhs;
        m_open_time = m_market->GetCurrentTime();
    }

    double GetCapacity() const {
        return m_balance_lhs + m_balance_rhs;
    }

    double GetOpportunityCost() const {

        double time_difference = m_market->GetCurrentTime() - m_open_time;

        auto capacity = GetCapacity();

        return m_market->BondValue(capacity, time_difference) - capacity;
    }

};
