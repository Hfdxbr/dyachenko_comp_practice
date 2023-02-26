#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <tuple>

template <class T, size_t N>
class Vector : public std::array<T, N> {
 public:
  Vector& operator+=(const Vector& other) {
    for (size_t i = 0; i < N; ++i)
      (*this)[i] += other[i];
    return *this;
  }

  Vector operator+(const Vector& other) const {
    Vector output(*this);
    output += other;
    return output;
  }

  Vector& operator-=(const Vector& other) {
    for (size_t i = 0; i < N; ++i)
      (*this)[i] -= other[i];
    return *this;
  }

  Vector operator-(const Vector& other) const {
    Vector output(*this);
    output -= other;
    return output;
  }

  Vector& operator*=(const Vector& other) {
    for (size_t i = 0; i < N; ++i)
      (*this)[i] *= other[i];
    return *this;
  }

  Vector operator*(const Vector& other) const {
    Vector output(*this);
    output *= other;
    return output;
  }

  Vector& operator*=(T k) {
    for (size_t i = 0; i < N; ++i)
      (*this)[i] *= k;
    return *this;
  }

  Vector operator*(T k) const {
    Vector output(*this);
    output *= k;
    return output;
  }

  auto as_tuple() const { return as_tuple(std::make_index_sequence<N>{}); }

  double norm1() const {
    return std::accumulate(this->begin(), this->end(), 0., [](auto norm, auto v) { return norm + std::abs(v); });
  }

  double norm2() const {
    return std::sqrt(std::accumulate(this->begin(), this->end(), 0., [](auto norm, auto v) { return norm + v * v; }));
  }

  static T dot(const Vector& v1, const Vector& v2) {
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
      sum += v1[i] * v2[i];
    return sum;
  }

 private:
  template <std::size_t... Is>
  auto as_tuple(std::index_sequence<Is...>) const {
    return std::tie((*this)[Is]...);
  }
};
