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

constexpr int N = 1000;
constexpr double tau = 1.0 / (N - 1);
constexpr double h = 1.0 / (N - 1);

struct Params {
  double alpha;
  double beta;
};

using Layer = std::array<double, N>;
using Solution = std::vector<Layer>;

double dist(const Layer& a, const Layer b) {
  double sum = 0;
  for (int i = 0; i < N; ++i) sum += std::pow(a[i] - b[i], 2);
  return std::sqrt(sum);
}

void tri_solve(const Layer& A, Layer& B, const Layer& C, Layer& F) {
  for (int i = 1; i < N; ++i) {
    double m = A[i - 1] / B[i - 1];
    B[i] -= m * C[i - 1];
    F[i] -= m * F[i - 1];
  }

  F[N - 1] /= B[N - 1];
  for (int i = N - 2; i >= 0; --i) {
    F[i] -= C[i] * F[i + 1];
    F[i] /= B[i];
  }
}

void fill(Layer& A, Layer& B, Layer& C, Layer& F, const Layer& u, double t, const Params& params) {}

void set_boundary(Layer& A, Layer& B, Layer& C, Layer& F, double t, const Params& params) {
  A.front() = C.front() = 0;
  B.front() = 1;
  F.front() = t > 0 ? t / (t + params.beta) : 0.0;
  A.back() = C.back() = 0;
  B.back() = 1;
  F.back() = 0;
}

Solution solve(const Params& params, double& error) {
  Solution sol = {Layer{}};
  sol.back().fill(0);
  Layer A, B, C, F;
  for (int i = 1; i < N; ++i) {
    double t = tau * i;
    Layer x;
    fill(A, B, C, F, sol.back(), t, params);
    set_boundary(A, B, C, F, t, params);
    tri_solve(A, B, C, F);
    do {
      x = std::move(F);
      fill(A, B, C, F, x, t, params);
      set_boundary(A, B, C, F, t, params);
      tri_solve(A, B, C, F);
    } while (dist(x, F) > eps);
    sol.push_back(F);
  }

  for (int i = 1; i < N - 1; ++i) {
    error += 0;
  }

  return sol;
}

void print_points(const Solution& sol, const std::string& filename) {
  std::ofstream ofs(filename);
  std::for_each(sol.begin(), sol.end(), [&ofs](const Layer& u) {
    ofs << u[0];
    for (int i = 1; i < N; ++i) ofs << ',' << u[i];
    ofs << std::endl;
  });
};

void execute(const Params& params) {
  bool stats_exists = std::ifstream("stats.csv").good();
  std::ofstream ofs_stats("stats.csv", std::ios_base::app);
  if (!stats_exists) ofs_stats << "\\(\\alpha\\),\\(\\beta\\),\\(N\\),\\(error\\)" << std::endl;
  double error = 0;
  Solution sol = solve(params, error);
  std::stringstream out_file;
  out_file << "points_" << to_string(params.alpha, 2, std::fixed) << "_"  //
           << to_string(params.beta, 2, std::fixed) << ".csv";
  print_points(sol, out_file.str());
  ofs_stats << to_string(params.alpha, 2, std::fixed) << ','  //
            << to_string(params.beta, 2, std::fixed) << ','   //
            << N << ',' << to_string(error, 3, std::scientific) << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "alpha and beta parameters required" << std::endl;
    return 1;
  }
  double alpha = std::stod(argv[1]);
  double beta = std::stod(argv[2]);
  execute({alpha, beta});

  return 0;
}