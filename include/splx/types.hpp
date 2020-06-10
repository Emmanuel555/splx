#ifndef SPLX_INTERNAL_TYPES_HPP
#define SPLX_INTERNAL_TYPES_HPP

#include <Eigen/Dense>

namespace splx {

template<typename T>
using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T, unsigned int DIM>
using VectorDIM = Eigen::Matrix<T, DIM, 1>;

template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T, unsigned int DIM>
using Hyperplane = Eigen::Hyperplane<T, DIM>;

using Index = Eigen::Index;

template<typename T, unsigned int DIM>
using AlignedBox = Eigen::AlignedBox<T, DIM>;

template<typename T>
struct Constraint {
    using _Row = Row<T>;

    Constraint(const _Row& c, T l, T u, bool sc = false, T sw = T(1)):
        coeff(c), lb(l), ub(u), soft_convertible(sc), soft_weight(sw) {

    }

    Constraint(const _Row&& c, T l, T u, bool sc = false, T sw = T(1)):
        coeff(std::move(c)), lb(l), ub(u), soft_convertible(sc), soft_weight(sw) {

    }

    _Row coeff;
    T lb;
    T ub;
    bool soft_convertible;
    T soft_weight;
};

}

#endif