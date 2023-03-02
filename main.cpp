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

double eps;
bool print_data = false;

using Vector2D = Vector<double, 2>;
using Point = Vector2D;
using Coeffs = Vector2D;
using TimePoint = std::pair<double, Point>;
using Solution = std::vector<TimePoint>;
using namespace std::string_literals;

double f1(const Point& p, double t) {
  const auto& [x1, x2] = p.as_tuple();
  return x2;
}

double f2(const Point& p, double t) {
  const auto& [x1, x2] = p.as_tuple();
  double t2 = t * t;
  return x1 / (1. + t2 * t2);
}

double phi1(const Solution& sol) {
  const auto& [x1, x2] = sol.back().second.as_tuple();
  return x1 - 1.;
}

double phi2(const Solution& sol) {
  const auto& [x1, x2] = sol.back().second.as_tuple();
  return x2 - 0.;
}

double norm(const Vector2D& v) { return v.norm1(); }

Coeffs calculate_k(Point p, double t, const Coeffs& k = {}) {
  p += k;
  return {f1(p, t), f2(p, t)};
}

std::pair<Point, Vector2D> gammas_with_error(const Point& p, double t, double h) {
  Coeffs k1 = calculate_k(p, t) * h;
  Coeffs k2 = calculate_k(p, t + h * 0.5, k1 * 0.5) * h;
  Coeffs k3 = calculate_k(p, t + h * 0.5, (k1 + k2) * 0.25) * h;
  Coeffs k4 = calculate_k(p, t + h, k3 * 2. - k2) * h;
  Coeffs k5 = calculate_k(p, t + h * (2. / 3.), (k1 * 7. + k2 * 10. + k4) * (1. / 27.)) * h;
  Coeffs k6 = calculate_k(p, t + h * 0.2, (k1 * 28. - k2 * 125. + k3 * 546. + k4 * 54. - k5 * 378.) * (1. / 625.)) * h;

  Point new_p = p + k1 * (1. / 24.) + k4 * (5. / 48.) + k5 * (27. / 56.) + k6 * (125. / 336.);
  Vector2D err = (k1 * 42. + k3 * 224. + k4 * 21. - k5 * 162. - k6 * 125.) * (-1. / 336.);
  return {new_p, err};
}

double update_step(const Vector2D& err, double h) {
  return std::clamp(std::pow(eps / (0.001 * eps + norm(err)), 1. / 6.), 0.1, 10.) * 0.95 * h;
}

Solution solve(double x10, double x20, double& error) {
  double t = 0;
  // double T = 1e7;
  double h = eps * 100.;
  Solution sol = {TimePoint{0, {x10, x20}}};
  while (f2(sol.back().second, t) > 0.1 * eps) {
    auto [p, err] = gammas_with_error(sol.back().second, t, h);
    if (double error = norm(err); error > eps || error < 0.1 * eps) {
      h = update_step(err, h);
      continue;
    }

    if (norm(sol.back().second - p) < eps)
      break;

    error += norm(err);
    t += h;
    sol.emplace_back(t, p);
  }
  return sol;
}

Solution shooting(double& x10, double& x20, double& error) {
  double dx10 = 0.001;
  double dx20 = 0.001;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-10.0, 10.0);

  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Solving " << std::endl;
  while (true) {
    std::cout << "\tparams: x10=" << x10 << ", x20=" << x20 << std::endl;
    Solution sol = solve(x10, x20, error);
    Vector2D phi{phi1(sol), phi2(sol)};

    if (phi.norm1() < eps) {
      std::cout << "Found!\n\n";
      reset_stream(std::cout);
      return sol;
    }

    double W11, W12, W21, W22;
    {
      double error = 0;
      Solution sol_ = solve(x10 + dx10, x20, error);
      Vector2D dphi = Vector2D{phi1(sol_), phi2(sol_)} - phi;
      W11 = dphi[0] / dx10;
      W21 = dphi[1] / dx10;
    }
    {
      double error = 0;
      Solution sol_ = solve(x10, x20 + dx20, error);
      Vector2D dphi = Vector2D{phi1(sol_), phi2(sol_)} - phi;
      W12 = dphi[0] / dx10;
      W22 = dphi[1] / dx10;
    }

    double detW = W11 * W22 - W12 * W21;
    if (detW < eps) {
      std::cout << "Singularity matrix can not be inversed, resettings params...\n";
      x10 = dis(gen);
      x20 = dis(gen);
      continue;
    }

    Vector2D new_params;
    new_params[0] = x10 - Vector2D::dot(Vector2D{W22, -W12}, phi) / detW;
    new_params[1] = x20 - Vector2D::dot(Vector2D{-W21, W11}, phi) / detW;

    std::tie(x10, x20) = new_params.as_tuple();
  }
}

void print_points(const Solution& sol, const std::string& filename) {
  std::ofstream ofs(filename);
  PrintHelper ph(ofs, "", "\n", ",");
  ph.print("t", "x1", "x2");
  std::for_each(sol.begin(), sol.end(), [&ph](const TimePoint& tp) {
    const auto& [t, p] = tp;
    const auto& [x1, x2] = p.as_tuple();
    ph.print(t, x1, x2);
  });
};

void execute(double x10, double x20) {
  bool stats_exists = std::ifstream("stats.csv").good();
  std::ofstream ofs_stats("stats.csv", std::ios_base::app);
  ofs_stats << std::setprecision(8) << std::fixed;
  PrintHelper ph(ofs_stats, "", "\n", ",");
  if (!stats_exists)
    ph.print("\\(\\varepsilon\\)", "n", "\\(x_1(0)\\)", "\\(x_2(0)\\)", "\\(x_1(\\infty)\\)", "\\(error\\)");

  double error = 0;
  Solution sol = shooting(x10, x20, error);

  if (print_data)
    print_points(sol, "data.csv");
  const auto [x1, x2] = sol.back().second.as_tuple();
  ph.print(to_string(eps, 1, std::scientific), sol.size(), x10, x20, x1, to_string(error, 3, std::scientific));
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "eps parameters required" << std::endl;
    return 1;
  }
  eps = std::stod(argv[1]);
  if (argc == 3 && argv[2] == "print"s)
    print_data = true;
  double x10 = 1.0;
  double x20 = -1.0;
  execute(x10, x20);

  return 0;
}