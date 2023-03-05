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
double beta;

double f1(const Points& u, double h) {
  double sum = 0;
  for (int i = 1; i < u.size() - 1; i++)
    sum += u[i] * h;
  sum += 0.5 * (u[0] * h + u.back() * h);
  return sum;
}

double f2(const Points& u, double h) {
  double sum = 0;
  for (int i = 1; i < u.size() - 1; i++)
    sum += std::pow(u[i], 4.) * h;
  sum += 0.5 * (std::pow(u[0], 4.) * h + std::pow(u.back(), 4.) * h);
  return sum - 50 * alpha * std::pow(u.back(), 4.);
}

double A(double h, double tau) { return -1. / (2. * h * h); }
double B(double h, double tau) { return 1. / (h * h) + 1. / tau; }
double B_hat(double h, double tau) { return -1. / (h * h) + 1. / tau; }
double C(double h, double tau) { return -1. / (2. * h * h); }

void fill_f(const Points& u, double h, double tau, Points& f) {
  for (int i = 1; i < f.size() - 1; ++i)
    f[i] = -A(h, tau) * u[i - 1] + B_hat(h, tau) * u[i] - C(h, tau) * u[i + 1] + std::pow(u[i], 4.);

  int n = u.size() - 1;
  double r = (4. * h * std::pow(u[n], 4.) - (3. * u[n] - 4. * u[n - 1] + u[n - 2])) / 3.;
  f[n - 1] += r * C(h, tau);
}

void fill_diags(double h, double tau, Points& As, Points& Bs, Points& Cs) {
  for (int i = 1; i < As.size() - 1; ++i) {
    As[i] = A(h, tau);
    Bs[i] = B(h, tau);
    Cs[i] = C(h, tau);
  }
  Bs[1] += As[1] * 4. / 3.;
  Cs[1] -= As[1] / 3.;
  As[1] = 0;

  int n = As.size() - 1;
  Bs[n] += Cs[n] * 4. / 3.;
  As[n] -= Cs[n] / 3.;
}

void solve_tri(const Points& As, Points Bs, const Points& Cs, Points& f, Points& u) {
  int n = f.size() - 1;
  for (int i = 2; i < n; ++i) {
    double m = As[i] / Bs[i - 1];
    Bs[i] = Bs[i] - m * Cs[i - 1];
    f[i] = f[i] - m * f[i - 1];
  }

  u[n - 1] = f[n - 1] / Bs[n - 1];

  for (int i = n - 2; i > 0; i--) {
    u[i] = (f[i] - Cs[i] * u[i + 1]) / Bs[i];
  }
}

void solve_next(double h, double tau, Points& u) {
  static Points As(u.size());
  static Points Bs(u.size());
  static Points Cs(u.size());
  static Points f(u.size());
  fill_diags(h, tau, As, Bs, Cs);
  fill_f(u, h, tau, f);

  int n = u.size() - 1;
  double r = (4. * h * std::pow(u[n], 4.) - (3. * u[n] - 4. * u[n - 1] + u[n - 2])) / 3.;

  solve_tri(As, Bs, Cs, f, u);
  u[0] = (4 * u[1] - u[2]) / 3.0;
  u[n] = r + u[n - 1] * 4. / 3. - u[n - 2] / 3.;
}

void process(int N, int M) {
  double h = 1.0 / (N - 1);
  double tau = 1.0 / (M - 1);
  Points u(N);
  auto e = "_a" + to_string(alpha, 1) + "_b" + to_string(beta, 1) + "_N" + std::to_string(N) + "_M" + std::to_string(M);
  std::ofstream ofs("fs" + e + ".csv");
  std::ofstream ou("u" + e + ".csv");
  ofs << "t,f1,f2\n";
  ou << "\\(t\\)";
  for (auto i : {"0.0", "0.2", "0.4", "0.6", "0.8", "1.0"})
    ou << ",\\(\\left.u\\right\\vert_{x=" << i << "}\\)";
  ou << "\n";
  for (int i = 0; i < N; ++i)
    u[i] = beta * std::pow(1. - std::pow(i * h, 2.), 2);

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
  beta = std::stod(argv[4]);
  process(N, M);
  return 0;
}