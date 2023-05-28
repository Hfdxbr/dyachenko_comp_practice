#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

constexpr double pi = M_PI;

constexpr double a_k(int k) {
  constexpr int fact[]{1, 1, 2, 6, 24};
  return (pi - k) * (pi - k) / fact[k];
}

double f_k(double x, int k) {
  constexpr double as[]{a_k(0), a_k(1), a_k(2), a_k(3)};
  return as[k] * std::pow(x, 2 * k + 1);
}

double g(double x, double alpha) {
  double s = 0.1 * alpha;
  for (int k = 0; k <= 3; ++k) s -= f_k(x, k);
  return s;
}

double phi1(double x, double alpha) {
  constexpr double a_0 = a_k(0);
  double s = 0.1 * alpha;
  for (int k = 1; k <= 3; ++k) s -= f_k(x, k);
  return s / a_0;
}

double phi2(double x, double alpha) {
  constexpr double ks[]{5e-1, 5e-2, 5e-3, 5e-4, 5e-5, 5e-6};
  int n = static_cast<int>(std::log10(alpha));
  return x + ks[n] * g(x, alpha);
}

double solve(double alpha, double eps, double x0) {
  auto func = (std::abs(alpha) <= 100) ? phi1 : phi2;
  double x = func(x0, alpha);
  while (std::abs(g(x, alpha)) > eps) {
    if (std::abs(x - x0) < std::numeric_limits<double>::epsilon() * std::max(x, x0)) break;
    x0 = x;
    x = func(x0, alpha);
  }
  return x;
}

int main() {
  std::array alphas{0.1, 1., 10., 100., 1000., 10000., 100000.};
  std::array epss{1e-7, 1e-9, 1e-11};
  std::ofstream ofs("data/stats.csv");
  double x0 = 0.0;
  ofs << "\\(\\alpha\\),\\(\\varepsilon\\),x,g(x)" << std::setprecision(13) << std::endl;
  for (auto alpha : alphas) {
    for (auto eps : epss) {
      double x = solve(alpha, eps, x0);
      ofs << alpha << ',' << eps << ',' << x << ',' << g(x, alpha) << std::endl;
    }
  }
  return 0;
}