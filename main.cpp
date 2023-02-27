#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "printer.h"
#include "vector.h"

constexpr double eps = 1e-11;

using Vector2D = Vector<double, 2>;
using Point = Vector2D;
using Coeffs = Vector2D;
using TimePoint = std::pair<double, Point>;
using Solution = std::vector<TimePoint>;

double f1(const Point& p, double t, double alpha) {
  const auto& [x1, x2] = p.as_tuple();
  return x2;
}

double f2(const Point& p, double t, double alpha) {
  const auto& [x1, x2] = p.as_tuple();
  return (1 + t * t) * x1 / alpha;
}

double phi(const Solution& sol, double alpha) {
  const auto& [x1, x2] = sol.back().second.as_tuple();
  return x2 - 1.;
}

double norm(const Vector2D& v) { return v.norm1(); }

Coeffs calculate_k(Point p, double t, double alpha, const Coeffs& k = {}) {
  p += k;
  return {f1(p, t, alpha), f2(p, t, alpha)};
}

std::pair<Point, Vector2D> gammas_with_error(const Point& p, double t, double alpha, double h) {
  Coeffs k1 = calculate_k(p, t, alpha) * h;
  Coeffs k2 = calculate_k(p, t + h * 0.5, alpha, k1 * 0.5) * h;
  Coeffs k3 = calculate_k(p, t + h * 0.5, alpha, (k1 + k2) * 0.25) * h;
  Coeffs k4 = calculate_k(p, t + h, alpha, k3 * 2. - k2) * h;
  Coeffs k5 = calculate_k(p, t + h * (2. / 3.), alpha, (k1 * 7. + k2 * 10. + k4) * (1. / 27.)) * h;
  Coeffs k6 =
    calculate_k(p, t + h * 0.2, alpha, (k1 * 28. - k2 * 125. + k3 * 546. + k4 * 54. - k5 * 378.) * (1. / 625.)) * h;

  Point new_p = p + k1 * (1. / 24.) + k4 * (5. / 48.) + k5 * (27. / 56.) + k6 * (125. / 336.);
  Vector2D err = (k1 * 42. + k3 * 224. + k4 * 21. - k5 * 162. - k6 * 125.) * (-1. / 336.);
  return {new_p, err};
}

double update_step(const Vector2D& err, double h) {
  return std::clamp(std::pow(eps / (0.001 * eps + norm(err)), 1. / 6.), 0.1, 10.) * 0.95 * h;
}

Solution solve(double x10, double alpha, double beta, double& error) {
  double t = 0;
  double T = beta;
  double h = eps * 100.;
  double max_allowed_h = T / 100;  // For smoother graphic
  Solution sol = {TimePoint{0, {x10, 0}}};
  while (t < T) {
    auto [p, err] = gammas_with_error(sol.back().second, t, alpha, h);
    if (double error = norm(err); error > eps || error < 0.1 * eps) {
      h = update_step(err, h);
      continue;
    }

    bool need_recalc = false;
    if (h > max_allowed_h) {
      h = max_allowed_h;
      need_recalc = true;
    }

    if (t + h >= T) {
      h = T - t;
      need_recalc = true;
    }

    if (need_recalc)
      std::tie(p, err) = gammas_with_error(sol.back().second, t, alpha, h);

    error += norm(err);
    t += h;
    sol.emplace_back(t, p);
  }
  return sol;
}

Solution shooting(double& x10, double alpha, double beta, double& error) {
  double dx10 = 0.001;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-10.0, 10.0);

  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Solving for alpha=" << to_string(alpha, 1) << ", beta=" << to_string(beta, 1) << std::endl;
  while (true) {
    std::cout << "\tparams: x10=" << x10 << std::endl;
    Solution sol = solve(x10, alpha, beta, error);
    double phi_ = phi(sol, alpha);

    if (std::abs(phi_) < eps) {
      std::cout << "Found!\n\n";
      reset_stream(std::cout);
      return sol;
    }

    double dphi;
    {
      double error = 0;
      Solution sol_ = solve(x10 + dx10, alpha, beta, error);
      dphi = (phi(sol_, alpha) - phi_) / dx10;
    }

    if (std::abs(dphi) < eps) {
      std::cout << "Division by zero occured, resettings params...\n";
      x10 = dis(gen);
      continue;
    }

    x10 = x10 - phi_ / dphi;
  }
}

void print_points(const Solution& sol, double alpha, const std::string& filename) {
  std::ofstream ofs(filename);
  PrintHelper ph(ofs, "", "\n", ",");
  ph.print("t", "x1", "x2");
  std::for_each(sol.begin(), sol.end(), [&ph, alpha](const TimePoint& tp) {
    const auto& [t, p] = tp;
    const auto& [x1, x2] = p.as_tuple();
    ph.print(t, x1, x2);
  });
};

void execute(double x10, double alpha, double beta) {
  bool stats_exists = std::ifstream("stats.csv").good();
  std::ofstream ofs_stats("stats.csv", std::ios_base::app);
  ofs_stats << std::setprecision(8) << std::fixed;
  PrintHelper ph(ofs_stats, "", "\n", ",");
  if (!stats_exists)
    ph.print("\\(\\alpha\\)", "\\(\\beta\\)", "\\(x_1(0)\\)", "\\(x_2(\\beta)\\)", "\\(error\\)");

  double error = 0;
  Solution sol = shooting(x10, alpha, beta, error);

  int ap = alpha > 0.1 ? 1 : 2;
  std::string plot_file = "a" + to_string(alpha, ap) + "_b" + to_string(beta, 1) + ".csv";
  print_points(sol, alpha, plot_file);
  const auto [x1, x2] = sol.back().second.as_tuple();
  ph.print(to_string(alpha, ap), to_string(beta, 1), x10, x2, to_string(error, 3, std::scientific));
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "alpha and beta parameters required" << std::endl;
    return 1;
  }
  double alpha = std::stod(argv[1]);
  double beta = std::stod(argv[2]);
  double x10 = argc == 4 ? std::stod(argv[3]) : 1.0;
  execute(x10, alpha, beta);

  return 0;
}