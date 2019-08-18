#include <cstdlib>
#include <cassert>
#include <vector>
#include <tuple>

template<typename T>
constexpr T ext_gcd(const T a, const T b, T* x, T* y) {
  if(b == static_cast<T>(0)) {
    *x = 1; *y = 0;
    return a;
  }
  const T d = ext_gcd(b, a % b, y, x);
  *y -= *x * (a / b);
  return d;
}

template<typename T>
constexpr std::vector<std::tuple<T, T>> prime_factorization(T n) {
  std::vector<std::tuple<T, T>> ret;

  T p = 2;
  while (n >= p * p) {
    T index = static_cast<T>(0);
    while (n % p == static_cast<T>(0)) {index++; n /= p;}
    if (index > static_cast<T>(0)) ret.emplace_back(p, index);
    p++;
  }
  if (n > static_cast<T>(1)) {
    ret.emplace_back(n, static_cast<T>(1));
  }
  return ret;
}

template <typename T>
T fast_pow(T a, T n) {
  T ret = static_cast<T>(1);
  while (n > static_cast<T>(0)) {
    if (n & static_cast<T>(1)) ret = ret * a;
    a *= a;
    n >>= static_cast<T>(1);
  }
  return ret;
}

template <typename T>
T mod(const T a, const T m) {
  T ret = a % m;
  return (ret >= 0) ? ret : ret + m;
}

int main(void) {
  int n;
  scanf("%d", &n);
  const auto ps = prime_factorization(n);
  for (const auto v : ps) {
    const int p = std::get<0>(v);
    const int e = std::get<1>(v);
    printf("%d^%d = %d ",p, e, fast_pow(p, e));
  }
  printf("\n\n");

  if (ps.size() < 2) {
    printf("Please select other n\n");
    return 0;
  }

  int n1 = fast_pow(std::get<0>(ps[0]), std::get<1>(ps[0]));
  assert (n % n1 == 0);
  int n2 = n / n1;

  printf("n1 = %d n2 = %d  n = %d\n", n1, n2, n1 * n2);

  int inv_n1,inv_n2;
  ext_gcd(n1, n2, &inv_n1, &inv_n2);
  inv_n1 = mod(inv_n1, n2);
  inv_n2 = mod(inv_n2, n1);

  printf ("n1^(-1)=%d n1xn1^(-1)=%d mod %d\n", inv_n1, mod(n1 * inv_n1, n2), n2);
  printf ("n2^(-1)=%d n2xn2^(-1)=%d mod %d\n", inv_n2, mod(n2 * inv_n2, n1), n1);

  assert(mod(n1 * inv_n1, n2) == 1);
  assert(mod(n2 * inv_n2, n1) == 1);

  return 0;
}