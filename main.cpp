#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "printer.h"
#include "vector.h"

using Points = std::vector<double>;
using TimePoint = std::pair<double, Points>;
using namespace std::string_literals;

double alpha;

double f1(const Points& u, double h) {
  double sum = 0;
  for (int i = 1; i < u.size() - 1; i++)
    sum += u[i] * (i * h) * h;
  sum += 0.5 * (u[0] * (0 * h) * h + u.back() * (1 * h) * h);
  return sum;
}

double f2(const Points& u, double h) {
  int n = u.size() - 1;
  return (3 * u[n] - 4 * u[n - 1] + u[n - 2]) / (2.0 * h);
}

double A(int i, double h, double tau) { return -(i - 0.5) / (2. * h * h * i); }

double B(int i, double h, double tau) { return 1.0 / tau + 1.0 / (h * h) + 0.5 * alpha; }

double B_hat(int i, double h, double tau) { return 1.0 / tau - 1.0 / (h * h) + 0.5 * alpha; }

double C(int i, double h, double tau) { return -(i + 0.5) / (2. * h * h * i); }

void fill_f(const Points& u_prev, double h, double tau, Points& f) {
  for (int i = 1; i < f.size() - 1; ++i)
    f[i] = -A(i, h, tau) * u_prev[i - 1] + B_hat(i, h, tau) * u_prev[i] - C(i, h, tau) * u_prev[i + 1] + 1;
  // u(1) = 0
  f.back() = 0;
}

void fill_diags(double h, double tau, Points& As, Points& Bs, Points& Cs) {
  for (int i = 1; i < As.size() - 1; ++i) {
    As[i] = A(i, h, tau);
    Bs[i] = B(i, h, tau);
    Cs[i] = C(i, h, tau);
  }
  Bs[1] += As[1] * 4. / 3.;
  Cs[1] -= As[1] / 3.;
  As[1] = 0;
  As.back() = 0.;
  Cs.back() = 0.;
  Bs.back() = 1.;
}

void solve_tri(const Points& As, Points Bs, const Points& Cs, Points& f, Points& u) {
  int n = f.size();
  for (int i = 2; i < n; ++i) {
    double m = As[i] / Bs[i - 1];
    Bs[i] = Bs[i] - m * Cs[i - 1];
    f[i] = f[i] - m * f[i - 1];
  }

  u[n - 1] = f[n - 1] / Bs[n - 1];

  for (int i = n - 2; i > 0; i--) {
    u[i] = (f[i] - Cs[i] * u[i + 1]) / Bs[i];
  }
  u[0] = (4 * u[1] - u[2]) / 3.0;
}

void solve_next(double h, double tau, Points& u) {
  static Points As(u.size());
  static Points Bs(u.size());
  static Points Cs(u.size());
  static Points f(u.size());
  fill_diags(h, tau, As, Bs, Cs);
  fill_f(u, h, tau, f);
  solve_tri(As, Bs, Cs, f, u);
}

void process(int N, int M) {
  double h = 1.0 / (N - 1);
  double tau = 1.0 / (M - 1);
  Points u(N);

  std::ofstream ofs("fs_a" + to_string(alpha, 1) + "_N" + std::to_string(N) + "_M" + std::to_string(M) + ".csv");
  std::ofstream ou("u_a" + to_string(alpha, 1) + "_N" + std::to_string(N) + "_M" + std::to_string(M) + ".csv");
  ofs << "t,f1,f2\n";
  ou << "\\(t\\),\\(\\left.u\\right\\vert_{x=0.0}\\),\\(\\left.u\\right\\vert_{x=0.2}\\),\\(\\left.u\\right\\vert_{x=0.4}\\),\\(\\left.u\\right\\vert_{x=0.6}\\),\\(\\left.u\\right\\vert_{x=0.8}\\),\\(\\left.u\\right\\vert_{x=1.0}\\)\n";
  for (int i = 0; i < N; ++i)
    u[i] = 0.5 * (1. - i * i * h * h);

  ofs << 0 << "," << f1(u, h) << "," << f2(u, h) << "\n";
  ou << 0;
  for (int i = 0; i < N; i += N / 5)
    ou << "," << u[i];
  ou << "\n";

  for (int j = 1; j < M; ++j) {
    solve_next(h, tau, u);

    ofs << j * tau << "," << f1(u, h) << "," << f2(u, h) << "\n";
    if (j % (M / 10) == 0) {
      ou << j * tau;
      for (int i = 0; i < N; i += N / 5)
        ou << "," << u[i];
      ou << "\n";
    }
  }
}

int main(int argc, char* argv[]) {
  int N = std::stoi(argv[1]);
  int M = std::stoi(argv[2]);
  alpha = std::stod(argv[3]);
  process(N, M);
  return 0;
}