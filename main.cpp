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
constexpr double T = 1.;
constexpr double max_allowed_h = T / 100;

using Vector2D = Vector<double, 2>;
using Vector4D = Vector<double, 4>;
using Point = Vector4D;
using Coeffs = Vector4D;
using TimePoint = std::pair<double, Point>;
using Solution = std::vector<TimePoint>;

double f1(const Point& p, double t, double alpha) {
  const auto& [x1, x2, p1, p2] = p.as_tuple();
  return x2;
}

double f2(const Point& p, double t, double alpha) {
  const auto& [x1, x2, p1, p2] = p.as_tuple();
  return p2;
}

double f3(const Point& p, double t, double alpha) {
  const auto& [x1, x2, p1, p2] = p.as_tuple();
  double denom = 2. + std::cos(alpha * x1);
  return -24. * alpha * x2 * std::sin(alpha * x1) / denom / denom;
}

double f4(const Point& p, double t, double alpha) {
  const auto& [x1, x2, p1, p2] = p.as_tuple();
  double denom = 2. + std::cos(alpha * x1);
  return -24. / denom - p1;
}

double L(const Point& p, double t, double alpha) {
  const auto& [x1, x2, p1, p2] = p.as_tuple();
  return p2 * p2 - 48. * x2 / (2. + std::cos(alpha * x1));
}

double phi1(const Solution& sol, double alpha) {
  const auto& [x1, x2, p1, p2] = sol.back().second.as_tuple();
  return x1 - 0.;
}

double phi2(const Solution& sol, double alpha) {
  const auto& [x1, x2, p1, p2] = sol.back().second.as_tuple();
  return p2 - 0.;
}

double J(const Solution& sol, double alpha, double p20) {
  auto sum = [L_prev = (p20 * p20), t_prev = 0., alpha](double s, const TimePoint& tp) mutable {
    const auto& [t_curr, p] = tp;
    double L_curr = L(p, t_curr, alpha);
    s += (L_curr + L_prev) * (t_curr - t_prev) * 0.5;
    t_prev = t_curr;
    L_prev = L_curr;
    return s;
  };

  return std::accumulate(std::next(sol.begin()), sol.end(), 0., sum);
}

double norm(const Vector4D& v) { return v.norm1(); }

Coeffs calculate_k(Point p, double t, double alpha, const Coeffs& k = {}) {
  p += k;
  return {f1(p, t, alpha), f2(p, t, alpha), f3(p, t, alpha), f4(p, t, alpha)};
}

std::pair<Point, Vector4D> gammas_with_error(const Point& p, double t, double alpha, double h) {
  Coeffs k1 = calculate_k(p, t, alpha) * h;
  Coeffs k2 = calculate_k(p, t + h * 0.5, alpha, k1 * 0.5) * h;
  Coeffs k3 = calculate_k(p, t + h * 0.5, alpha, (k1 + k2) * 0.25) * h;
  Coeffs k4 = calculate_k(p, t + h, alpha, k3 * 2. - k2) * h;
  Coeffs k5 = calculate_k(p, t + h * (2. / 3.), alpha, (k1 * 7. + k2 * 10. + k4) * (1. / 27.)) * h;
  Coeffs k6 =
    calculate_k(p, t + h * 0.2, alpha, (k1 * 28. - k2 * 125. + k3 * 546. + k4 * 54. - k5 * 378.) * (1. / 625.)) * h;

  Point new_p = p + k1 * (1. / 24.) + k4 * (5. / 48.) + k5 * (27. / 56.) + k6 * (125. / 336.);
  Vector4D err = (k1 * 42. + k3 * 224. + k4 * 21. - k5 * 162. - k6 * 125.) * (-1. / 336.);
  return {new_p, err};
}

constexpr double one_six = 1. / 6.;

double update_step(const Vector4D& err, double h) {
  return std::clamp(std::pow(eps / (0.001 * eps + norm(err)), one_six), 0.1, 10.) * 0.95 * h;
}

Solution solve(double x10, double p20, double alpha, double& error) {
  double t = 0;
  double h = eps * 100.;
  Solution sol = {TimePoint{0, {x10, 0, 0, p20}}};
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

Solution shooting(double& x10, double& p20, double alpha, double& error) {
  double dx10 = 0.001;
  double dp20 = 0.001;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-10.0, 10.0);

  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Solving for alpha=" << to_string(alpha, 1) << std::endl;
  while (true) {
    std::cout << "\tparams: x10=" << x10 << ", p20=" << p20 << std::endl;
    Solution sol = solve(x10, p20, alpha, error);
    Vector2D phi{phi1(sol, alpha), phi2(sol, alpha)};

    if (phi.norm1() < eps) {
      std::cout << "Found!\n\n";
      reset_stream(std::cout);
      return sol;
    }

    double W11, W12, W21, W22;
    {
      double error = 0;
      Solution sol_ = solve(x10 + dx10, p20, alpha, error);
      Vector2D dphi = Vector2D{phi1(sol_, alpha), phi2(sol_, alpha)} - phi;
      W11 = dphi[0] / dx10;
      W21 = dphi[1] / dx10;
    }
    {
      double error = 0;
      Solution sol_ = solve(x10, p20 + dp20, alpha, error);
      Vector2D dphi = Vector2D{phi1(sol_, alpha), phi2(sol_, alpha)} - phi;
      W12 = dphi[0] / dx10;
      W22 = dphi[1] / dx10;
    }

    double detW = W11 * W22 - W12 * W21;
    if (detW < eps) {
      std::cout << "Singularity matrix can not be inversed, resettings params...\n";
      x10 = dis(gen);
      p20 = dis(gen);
      continue;
    }

    Vector2D new_params;
    new_params[0] = x10 - Vector2D::dot(Vector2D{W22, -W12}, phi) / detW;
    new_params[1] = p20 - Vector2D::dot(Vector2D{-W21, W11}, phi) / detW;

    std::tie(x10, p20) = new_params.as_tuple();
  }
}

void print_points(const Solution& sol, double alpha, const std::string& filename) {
  std::ofstream ofs(filename);
  PrintHelper ph(ofs, "", "\n", ",");
  ph.print("t", "x1", "x2", "p1", "p2", "L");
  std::for_each(sol.begin(), sol.end(), [&ph, alpha](const TimePoint& tp) {
    const auto& [t, p] = tp;
    const auto& [x1, x2, p1, p2] = p.as_tuple();
    ph.print(t, x1, x2, p1, p2, L(p, t, alpha));
  });
};

void execute(double x10, double p20, std::vector<double> alphas) {
  std::sort(alphas.begin(), alphas.end());

  std::ofstream ofs_stats("stats.csv");
  ofs_stats << std::setprecision(8) << std::fixed;
  PrintHelper ph(ofs_stats, "", "\n", ",");
  ph.print("\\(\\alpha\\)", "\\(x_1(0)\\)", "\\(p_2(0)\\)", "\\(x_2(1)\\)", "\\(p_1(1)\\)", "\\(J\\)", "\\(error\\)");

  std::for_each(alphas.begin(), alphas.end(), [&x10, &p20, &ph](double alpha) {
    double error = 0;
    Solution sol = shooting(x10, p20, alpha, error);

    std::string plot_file = "points_alpha" + to_string(alpha, 1) + ".csv";
    print_points(sol, alpha, plot_file);
    const auto [x1, x2, p1, p2] = sol.back().second.as_tuple();
    ph.print(to_string(alpha, 1), x10, p20, x2, p1, J(sol, alpha, p20), to_string(error, 3, std::scientific));
  });
}

int main() {
  std::vector<double> alphas{0.0, 0.1, 1.0, 5.1};
  double x10 = 1., p20 = -1.;
  execute(x10, p20, alphas);

  return 0;
}