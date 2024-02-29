#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

double eps;
bool print_data = false;

using Modifier = decltype(std::defaultfloat);

template <class T>
std::string to_string(const T& value, int precision = 6, Modifier modifier = std::fixed) {
  std::ostringstream out;
  out << modifier << std::setprecision(precision) << value;
  return out.str();
}

using Vector2D = std::array<double, 2>;
using Vector3D = std::array<double, 3>;
using Solution = std::vector<Vector3D>;


double f1(const Vector2D& p, double t, double lambda) {
  const auto& [x1, x2] = p;
  return x2;
}

double f2(const Vector2D& p, double t, double lambda) {
  const auto& [x1, x2] = p;
  return 1. - lambda * std::pow(t, 2) - std::pow(x1, 2);
}

double phi(const Solution& sol, double lambda) {
  const auto& [t, x1, x2] = sol.back();
  return x1 - 0.;
}

double norm(const Vector2D& v) {
  const auto& [v1, v2] = v;
  return std::sqrt(v1 * v1 + v2 * v2);
}

Vector2D tmp;
template <class F>  // F = update_tmp_ith_pos(int i);
void calculate_k_inplace(Vector2D& k, double t, double lambda, double h, F&& func) {
  for (int i = 0; i < tmp.size(); ++i)
    tmp[i] = func(i);
  k = {f1(tmp, t, lambda) * h, f2(tmp, t, lambda) * h};
}

Vector2D k1, k2, k3, k4, k5, k6;
Vector2D p, err;
void gammas_with_error_inplace(const Vector3D& tp, double lambda, double h){
  auto t = tp[0];
  calculate_k_inplace(k1, t, lambda, h, [&tp](int i) { return tp[i+1]; });
  calculate_k_inplace(k2, t + h * 0.5, lambda, h, [&tp](int i) { return tp[i+1] + k1[i] * 0.5; });
  calculate_k_inplace(k3, t + h * 0.5, lambda, h, [&tp](int i) { return tp[i+1] + (k1[i] + k2[i]) * 0.25; });
  calculate_k_inplace(k4, t + h, lambda, h, [&tp](int i) { return tp[i+1] + k3[i] * 2 - k2[i]; });
  calculate_k_inplace(k5, t + h * (2. / 3.), lambda, h, [&tp](int i) { return tp[i+1] + (k1[i] * 7 + k2[i] * 10 + k4[i]) / 27; });
  calculate_k_inplace(k6, t + h * 0.2, lambda, h, [&tp](int i) { return tp[i+1] + (k1[i] * 28 - k2[i] * 125 + k3[i] * 546 + k4[i] * 54 - k5[i] * 378) / 625; });

  for (int i = 0; i < p.size(); ++i){
    p[i] = tp[i+1] + k1[i] / 24 + k4[i] * (5. / 48) + k5[i] * (27. / 56) + k6[i] * (125. / 336);
    err[i] = -(k1[i] * 42 + k3[i] * 224 + k4[i] * 21 - k5[i] * 162 - k6[i] * 125) / 336;
  }
}

double update_step(const Vector2D& err, double h) {
  return std::clamp(std::pow(eps / (0.001 * eps + norm(err)), 1. / 6.), 0.1, 10.) * 0.95 * h;
}

Solution solve(double lambda, double& error) {
  double t = 0.0;
  double T = 1.0;
  double h = eps * 100.;
  Solution sol = {Vector3D{0, 0, 1}};
  while (t < T) {
    gammas_with_error_inplace(sol.back(), lambda, h);
    if (double error = norm(err); error > eps || error < 0.1 * eps) {
      h = update_step(err, h);
      continue;
    }

    bool need_recalc = false;
    if (t + h >= T) {
      h = T - t;
      need_recalc = true;
    }

    if (need_recalc)
      gammas_with_error_inplace(sol.back(), lambda, h);

    error += norm(err);
    t += h;
    sol.push_back({t, p[0], p[1]});
  }
  return sol;
}

Solution shooting(double& lambda, double& error) {
  double delta = 0.001;

  while (true) {
    std::cout << "Solving for lambda=" << to_string(lambda, 10) << std::endl;
    Solution sol = solve(lambda, error);
    double phi_ = phi(sol, lambda);

    if (std::abs(phi_) < eps) {
      std::cout << "Found!" << std::endl;
      return sol;
    }

    double dphi;
    {
      double error = 0;
      Solution sol_ = solve(lambda + delta, error);
      dphi = (phi(sol_, lambda) - phi_) / delta;
    }

    if (std::abs(dphi) < eps)
      throw std::runtime_error("Division by zero occured, try another start point.");

    lambda = lambda - phi_ / dphi;
  }
}

void print_points(const Solution& sol, double lambda, const std::string& filename) {
  std::ofstream ofs(filename);
  ofs << "t,x1,x2" << std::endl;
  std::for_each(sol.begin(), sol.end(), [&ofs](const Vector3D& tp) {
    const auto& [t, x1, x2] = tp;
    ofs << t << ',' << x1 << ',' << x2 << std::endl;
  });
};

void execute(double lambda) {
  bool stats_exists = std::ifstream("stats.csv").good();
  std::ofstream ofs_stats("stats.csv", std::ios_base::app);
  if (!stats_exists)
    ofs_stats << "\\(\\varepsilon\\),\\(\\lambda\\),n,\\(x_1(1)\\),\\(error\\)" << std::endl;
  double error = 0;
  Solution sol = shooting(lambda, error);

  if (print_data)
    print_points(sol, lambda, "points.csv");
  const auto [t, x1, x2] = sol.back();
  ofs_stats << to_string(eps, 1, std::scientific) << ','  //
            << to_string(lambda, 10) << ','               //
            << sol.size() << ','                          //
            << x1 << ','                                  //
            << to_string(error, 3, std::scientific) << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "lambda and eps parameters required" << std::endl;
    return 1;
  }
  double lambda = std::stod(argv[1]);
  eps = std::stod(argv[2]);
  if (argc == 4 && argv[3] == std::string("print"))
    print_data = true;
  execute(lambda);

  return 0;
}