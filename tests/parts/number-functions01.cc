#include <cassert>
#include <cstdlib>
#include <tuple>
#include <vector>

template <typename T>
constexpr std::tuple<T, T> ext_gcd(const T a, const T b) {
  constexpr auto calc_xy = [](const T a, const T b, const std::tuple<T, T> xy) {
    const auto [x_dash, y_dash] = xy;
    return std::make_tuple(y_dash, x_dash - (a / b) * y_dash);
  };

  return (b == static_cast<T>(0)) ? std::make_tuple(1, 0)
                                  : calc_xy(a, b, ext_gcd(b, a % b));
}

template <typename T>
constexpr std::vector<std::tuple<T, T>> prime_factorization(T n) {
  std::vector<std::tuple<T, T>> ret;

  T p = 2;
  while (n >= p * p) {
    T index = static_cast<T>(0);
    while (n % p == static_cast<T>(0)) {
      index++;
      n /= p;
    }
    if (index > static_cast<T>(0)) ret.emplace_back(p, index);
    p++;
  }
  if (n > static_cast<T>(1)) {
    ret.emplace_back(n, static_cast<T>(1));
  }
  return ret;
}

template <typename T>
constexpr T fast_pow(T a, T n) {
  T ret = static_cast<T>(1);
  while (n > static_cast<T>(0)) {
    if (n & static_cast<T>(1)) ret = ret * a;
    a *= a;
    n >>= static_cast<T>(1);
  }
  return ret;
}

template <typename T>
constexpr T fast_mod_pow(T a, T n, const T m) {
  T ret = static_cast<T>(1);
  while (n > static_cast<T>(0)) {
    if (n & static_cast<T>(1)) ret = ret * a % m;
    a = a * a % m;
    n >>= static_cast<T>(1);
  }
  return ret;
}
template <typename T>
constexpr T mod(const T a, const T m) {
  T ret = a % m;
  return (ret >= 0) ? ret : ret + m;
}

template <typename T>
constexpr T euler_phi(const T n) {
  const std::vector<std::tuple<T, T>> factors = prime_factorization(n);
  T ret                                       = n;
  for (const auto [p, e] : factors) {
    ret *= (p - 1);
    ret /= p;
  }
  return ret;
}

template <typename T>
constexpr T find_one_of_primitive_root(const T p) {
  const std::vector<std::tuple<T, T>> factors = prime_factorization(p - 1);
  const T phi                                 = p - static_cast<T>(1);
  T ret                                       = static_cast<T>(2);
  // TODO : ret < p?
  for (; ret <= p; ret++) {
    bool ok = true;
    for (const auto [q, e] : factors) {
      if (fast_mod_pow(ret, phi / q, p) == 1) {
        ok = false;
        break;
      }
    }
    if (ok) break;
  }
  if (ret >= p) ret = static_cast<T>(-1);
  return ret;
}

int main(void) {
  int n;
  scanf("%d", &n);
  const auto ps = prime_factorization(n);
  for (const auto [p, e] : ps) {
    printf("%d^%d = %d ", p, e, fast_pow(p, e));
  }
  printf("\n\n");

  if (ps.size() < 2) {
    printf("Please select other n\n");
    return 0;
  }

  const int n1 = fast_pow(std::get<0>(ps[1]), std::get<1>(ps[1]));
  assert(n % n1 == 0);
  const int n2 = n / n1;

  printf("n1 = %d n2 = %d  n = %d\n", n1, n2, n1 * n2);

  auto [inv_n1, inv_n2] = ext_gcd(n1, n2);
  inv_n1                = mod(inv_n1, n2);
  inv_n2                = mod(inv_n2, n1);

  printf("n1^(-1)=%d n1xn1^(-1)=%d mod %d\n", inv_n1, mod(n1 * inv_n1, n2), n2);
  printf("n2^(-1)=%d n2xn2^(-1)=%d mod %d\n", inv_n2, mod(n2 * inv_n2, n1), n1);

  assert(mod(n1 * inv_n1, n2) == 1);
  assert(mod(n2 * inv_n2, n1) == 1);

  const int g1 = find_one_of_primitive_root(n1);

  printf("one of a primitive root is %d in mod %d\n", g1, n1);
  const int phi_n1 = euler_phi(n1);
  for (int i = 1; i < phi_n1; i++) {
    assert(fast_mod_pow(g1, i, n1) != 1 || i == n1 - 1);
  }
  printf("%d passed the primitive root test in mod %d\n", g1, n1);

  return 0;
}