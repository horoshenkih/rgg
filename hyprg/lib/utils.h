//
// Created by serkh on 11/12/16.
//

#ifndef HYPRG_UTILS_H
#define HYPRG_UTILS_H

#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <ctime>

#include "graph.h"

Graph read_graph_from_file(const char *filename);

//int random_i(int i) { return std::rand() % i; }
template <typename T> void shuffle_vector(vector<T> &vec) {
    std::random_shuffle(vec.begin(), vec.end());//, random_i);
}

template <class TRandomIterator, class OutputIterator>
OutputIterator n_random_elements(TRandomIterator begin, TRandomIterator end, OutputIterator out, unsigned n) {
    typedef typename std::iterator_traits<TRandomIterator>::value_type T;
    vector<T> copy_vec(begin, end);
    shuffle_vector(copy_vec);
    return std::copy(copy_vec.begin(), copy_vec.begin()+n, out);
}
#endif //HYPRG_UTILS_H
