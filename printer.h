#include <iomanip>
#include <iostream>

const auto default_precision = std::cout.precision();

class PrintHelper {
 public:
  PrintHelper(std::ostream& os, const std::string& prefix = "", const std::string& postfix = "",
              const std::string& sep = ",")
    : stream(os), prefix(prefix), sep(sep), postfix(postfix) {}

  PrintHelper& setSeparator(const std::string& sep) {
    this->sep = sep;
    return *this;
  }

  PrintHelper& setPrefix(const std::string& prefix) {
    this->prefix = prefix;
    return *this;
  }

  PrintHelper& setPostfix(const std::string& postfix) {
    this->postfix = postfix;
    return *this;
  }

  template <class F, typename... R>
  PrintHelper& print(F&& f, R&&... r) {
    stream << prefix << f;
    ((stream << sep << r), ...);
    stream << postfix;
    return *this;
  }

 private:
  std::ostream& stream;
  std::string prefix;
  std::string sep;
  std::string postfix;
};

void reset_stream(std::ostream& os) { os << std::setprecision(default_precision) << std::defaultfloat; }

using Modifier = decltype(std::defaultfloat);

template <class T>
std::string to_string(const T& value, int precision = 6, Modifier modifier = std::fixed) {
  std::ostringstream out;
  out << modifier << std::setprecision(precision) << value;
  return out.str();
}