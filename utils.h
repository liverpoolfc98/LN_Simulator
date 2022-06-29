#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <numeric>
#include <functional>
#include <iterator>

template <class Iterator>
class IteratorRange {
public:
    IteratorRange(Iterator begin, Iterator end) : begin_(begin), end_(end) {}

    Iterator begin() const { return begin_; }
    Iterator end() const { return end_; }

private:
    Iterator begin_, end_;
};

template <class Iterator, class Predicate>
class FilterIterator {
public:
    typedef typename std::iterator_traits<Iterator>::reference reference;
    FilterIterator(Iterator begin, Iterator end, Predicate predicate)
        : begin_(begin), end_(end), predicate_(predicate) {}
    Iterator Begin() {
        begin_ = std::find_if(begin_, end_, predicate_);
        return begin_;
    }
    FilterIterator& operator++() {
        if (begin_ == end_) {
            return *this;
        }
        ++begin_ = std::find_if(begin_, end_, predicate_);
        return *this;
    }
    reference operator*() const {
        return *begin_;
    }
    bool operator!=(FilterIterator elem) {
        return begin_ != elem.Begin();
    }
private:
    Iterator begin_, end_;
    Predicate predicate_;
};


template<class Predicate>
int BinSearch(int begin, int end, Predicate predicate) {
    // finds least iterator from [begin, end) such that predicate(iterator) = true
    // if there is no such iterator, returns end
    int left = begin, right = end;
    while (right > left) {
        int middle = left + (right - left) / 2;
        if (predicate(middle)) {
            right = middle;
        } else {
            left = middle + 1;
        }
    }
    return left;
}
