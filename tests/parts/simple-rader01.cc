#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <tuple>
#include <vector>

template <typename T>
constexpr std::vector<T> simple_dft(T* a, const uint32_t n,
                                    const T sign = -1.0) {
  std::vector<T> ret(n * 2);
  for (size_t k = 0; k < n; k++) {
    ret[k * 2] = ret[k * 2 + 1] = 0.0;
    for (size_t j = 0; j < n; j++) {
      const T theta   = sign * static_cast<T>(2.0 * M_PI) * k * j / n;
      const T omega_r = std::cos(theta);
      const T omega_i = std::sin(theta);
      ret[k * 2] += std::fma(a[j * 2], omega_r, -a[j * 2 + 1] * omega_i);
      ret[k * 2 + 1] += std::fma(a[j * 2], omega_i, a[j * 2 + 1] * omega_r);
    }
  }
  return ret;
}

template <typename T>
constexpr std::tuple<T, T> ext_gcd(const T a, const T b) {
  constexpr auto calc_xy = [](const T a, const T b, const std::tuple<T, T> xy) {
    const auto [x_dash, y_dash] = xy;
    return std::make_tuple(y_dash, x_dash - (a / b) * y_dash);
  };

  return (b == static_cast<T>(0)) ? std::make_tuple<T, T>(1, 0)
                                  : calc_xy(a, b, ext_gcd(b, a % b));
}

template <typename T>
bool is_prime(T n) {
  bool ret = true;
  for (T i = 2; i * i <= n; i++) {
    if (n % i == 0) {
      ret = false;
      break;
    }
  }
  return ret;
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
constexpr T get_inverse(const T a, const T n) {
  // if(!is_prime(n)) return 0;
  const auto [x, y] = ext_gcd(a, n);
  (void)y;
  return mod(x, n);
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

template <typename T>
void cooley_turkey_fft_split_radix_frequency_inplace(T* const a,
                                                     const uint32_t n,
                                                     const T sign = -1.0) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  for (uint32_t m = n; m > 2; m >>= 1) {
    const uint32_t m4 = m / 4;

    const T theta = sign * 2.0 * M_PI / static_cast<T>(m);
    for (int k = m; k <= n; k <<= 2) {
      for (int i = k - m; i < n; i += 2 * k) {
        for (int j = 0; j < m4; j++) {
          const uint32_t id04 = 2 * (j + i);
          const uint32_t id14 = id04 + 2 * m4;
          const uint32_t id24 = id14 + 2 * m4;
          const uint32_t id34 = id24 + 2 * m4;

          const T qr = a[id04] - a[id24];
          const T qi = a[id04 + 1] - a[id24 + 1];

          a[id04] += a[id24];
          a[id04 + 1] += a[id24 + 1];

          const T sr = sign * (-a[id14 + 1] + a[id34 + 1]);
          const T si = sign * (a[id14] - a[id34]);

          a[id14] += a[id34];
          a[id14 + 1] += a[id34 + 1];

          const T omega1r = std::cos(theta * j);
          const T omega1i = std::sin(theta * j);
          const T omega2r = std::fma(static_cast<T>(2.0) * omega1r, omega1r,
                                     -static_cast<T>(1.0));
          const T omega2i = static_cast<T>(2.0) * omega1r * omega1i;
          const T omega3r = std::fma(omega1r, omega2r, -omega1i * omega2i);
          const T omega3i = std::fma(omega1r, omega2i, omega1i * omega2r);

          const T tmp_r2 = qr + sr;
          const T tmp_i2 = qi + si;
          a[id24]        = std::fma(tmp_r2, omega1r, -tmp_i2 * omega1i);
          a[id24 + 1]    = std::fma(tmp_r2, omega1i, tmp_i2 * omega1r);

          const T tmp_r3 = qr - sr;
          const T tmp_i3 = qi - si;
          a[id34]        = std::fma(tmp_r3, omega3r, -tmp_i3 * omega3i);
          a[id34 + 1]    = std::fma(tmp_r3, omega3i, tmp_i3 * omega3r);
        }
      }
    }
  }

  for (uint32_t k = 2; k <= n; k <<= 2) {
    for (uint32_t j = k - 2; j < n; j += 2 * k) {
      const uint32_t id0 = 2 * j;
      const uint32_t id1 = 2 * j + 2;
      const T tmp_r      = a[id0] - a[id1];
      const T tmp_i      = a[id0 + 1] - a[id1 + 1];
      a[id0] += a[id1];
      a[id0 + 1] += a[id1 + 1];
      a[id1]     = tmp_r;
      a[id1 + 1] = tmp_i;
    }
  }
}

template <typename T>
void cooley_turkey_fft_split_radix_time_inplace(T* const a, const uint32_t n,
                                                const T sign = -1.0) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  for (uint32_t k = 2; k <= n; k <<= 2) {
    for (uint32_t j = k - 2; j < n; j += 2 * k) {
      const uint32_t id0 = 2 * j;
      const uint32_t id1 = 2 * j + 2;
      const T tmp_r      = a[id0] - a[id1];
      const T tmp_i      = a[id0 + 1] - a[id1 + 1];
      a[id0] += a[id1];
      a[id0 + 1] += a[id1 + 1];
      a[id1]     = tmp_r;
      a[id1 + 1] = tmp_i;
    }
  }

  for (uint32_t m = 4; m <= n; m <<= 1) {
    const uint32_t m4 = m / 4;

    const T theta = sign * 2.0 * M_PI / static_cast<T>(m);
    for (int j = m; j <= n; j <<= 2) {
      for (int i = j - m; i < n; i += 2 * j) {
        for (int k = 0; k < m4; k++) {
          const uint32_t id04 = 2 * (k + i);
          const uint32_t id14 = id04 + 2 * m4;
          const uint32_t id24 = id14 + 2 * m4;
          const uint32_t id34 = id24 + 2 * m4;

          const T omega1r = std::cos(theta * k);
          const T omega1i = std::sin(theta * k);
          const T omega2r = std::fma(static_cast<T>(2.0) * omega1r, omega1r,
                                     -static_cast<T>(1.0));
          const T omega2i = static_cast<T>(2.0) * omega1r * omega1i;
          const T omega3r = std::fma(omega1r, omega2r, -omega1i * omega2i);
          const T omega3i = std::fma(omega1r, omega2i, omega1i * omega2r);

          const T rr = std::fma(a[id24], omega1r, -a[id24 + 1] * omega1i);
          const T ri = std::fma(a[id24], omega1i, a[id24 + 1] * omega1r);
          const T sr = std::fma(a[id34], omega3r, -a[id34 + 1] * omega3i);
          const T si = std::fma(a[id34], omega3i, a[id34 + 1] * omega3r);

          const T xr = a[id04];
          const T xi = a[id04 + 1];
          const T yr = a[id14];
          const T yi = a[id14 + 1];
          const T zr = rr + sr;
          const T zi = ri + si;
          const T wr = sign * (-ri + si);
          const T wi = sign * (rr - sr);

          a[id04]     = xr + zr;
          a[id04 + 1] = xi + zi;
          a[id14]     = yr + wr;
          a[id14 + 1] = yi + wi;
          a[id24]     = xr - zr;
          a[id24 + 1] = xi - zi;
          a[id34]     = yr - wr;
          a[id34 + 1] = yi - wi;
        }
      }
    }
  }
}

template <typename T>
std::vector<T> cyclic_convolution_complex(const std::vector<T>& a,
                                          const std::vector<T>& b) {
  assert(a.size() == b.size());

  const size_t n       = a.size() / 2;
  const bool is_power2 = (!(n & (n - 1)));
  uint32_t N           = 2;
  while ((is_power2 ? 1 : 2) * n > N) N *= 2;

  std::vector<T> ac(2 * N);
  std::vector<T> bc(2 * N);
  for (size_t i = 0; i < N; i++) {
    ac[2 * i]     = (2 * i < a.size()) ? a[2 * i] : static_cast<T>(0.0);
    ac[2 * i + 1] = (2 * i + 1 < a.size()) ? a[2 * i + 1] : static_cast<T>(0.0);
    bc[2 * i]     = (2 * i < b.size()) ? b[2 * i] : static_cast<T>(0.0);
    bc[2 * i + 1] = (2 * i + 1 < b.size()) ? b[2 * i + 1] : static_cast<T>(0.0);
  }
  if (!is_power2) {
    for (size_t i = 0; i < n - 1; i++) {
      bc[2 * (N - i - 1)]     = b[2 * (n - i - 1)];
      bc[2 * (N - i - 1) + 1] = b[2 * (n - i - 1) + 1];
    }
  }
  cooley_turkey_fft_split_radix_frequency_inplace(ac.data(), N);
  cooley_turkey_fft_split_radix_frequency_inplace(bc.data(), N);
  for (size_t i = 0; i < N; i++) {
    const T i2    = 2 * i;
    const T tmp_r = ac[i2];
    const T tmp_i = ac[i2 + 1];
    ac[i2]        = tmp_r * bc[i2] - tmp_i * bc[i2 + 1];
    ac[i2 + 1]    = tmp_r * bc[i2 + 1] + tmp_i * bc[i2];
  }
  cooley_turkey_fft_split_radix_time_inplace(ac.data(), N, 1.0);
  std::vector<T> ret(2 * N);
  for (size_t i = 0; i < 2 * n; i++) {
    ret[i] = ac[i] / static_cast<T>(N);
  }
  return ret;
}

template <typename T>
void rader_fft(T* const a, const uint32_t n, const T sign = -1.0) {
  if (!is_prime(n)) return;

  const T a0r = a[0];
  const T a0i = a[1];

  for (uint32_t j = 1; j < n; j++) {
    a[0] += a[2 * j];
    a[1] += a[2 * j + 1];
  }

  const uint32_t g     = find_one_of_primitive_root(n);
  const uint32_t inv_g = get_inverse<int32_t>(g, n);

  std::vector<T> x((n - 1) * 2);
  std::vector<T> y((n - 1) * 2);
  uint32_t gq = 1;
  for (uint32_t q = 0; q < n - 1; q++) {
    x[2 * q]     = a[2 * gq];
    x[2 * q + 1] = a[2 * gq + 1];
    gq           = gq * g % n;
  }
  uint32_t inv_gq = 1;
  for (uint32_t q = 0; q < n - 1; q++) {
    const T theta = sign * static_cast<T>(2.0 * M_PI) * inv_gq / n;
    y[2 * q]      = cos(theta);
    y[2 * q + 1]  = sin(theta);
    inv_gq        = inv_gq * inv_g % n;
  }
  std::vector<T> convo = cyclic_convolution_complex(x, y);
  for (uint32_t p = 0; p < n - 1; p++) {
    convo[2 * p] += a0r;
    convo[2 * p + 1] += a0i;
  }

  uint32_t inv_gp = 1;
  for (uint32_t p = 0; p < n - 1; p++) {
    a[2 * inv_gp]     = convo[2 * p];
    a[2 * inv_gp + 1] = convo[2 * p + 1];
    inv_gp            = inv_gp * inv_g % n;
  }
}

template <typename T>
std::vector<T> convolution(const std::vector<T>& a, const std::vector<T>& b) {
  const size_t n = 2 * std::max(a.size(), b.size());
  int N          = n;
  while (!is_prime(N)) N++;
  std::vector<T> ac(2 * N);
  std::vector<T> bc(2 * N);
  for (size_t i = 0; i < N; i++) {
    ac[2 * i]     = (i < a.size()) ? a[i] : static_cast<T>(0.0);
    ac[2 * i + 1] = static_cast<T>(0.0);
    bc[2 * i]     = (i < b.size()) ? b[i] : static_cast<T>(0.0);
    bc[2 * i + 1] = static_cast<T>(0.0);
  }

  rader_fft(ac.data(), N);
  rader_fft(bc.data(), N);
  for (size_t i = 0; i < N; i++) {
    const T i2    = 2 * i;
    const T tmp_r = ac[i2];
    const T tmp_i = ac[i2 + 1];
    ac[i2]        = tmp_r * bc[i2] - tmp_i * bc[i2 + 1];
    ac[i2 + 1]    = tmp_r * bc[i2 + 1] + tmp_i * bc[i2];
  }
  rader_fft(ac.data(), N, 1.0);
  std::vector<T> ret(n);
  for (size_t i = 0; i < n; i++) {
    ret[i] = (ac[i * 2] / static_cast<T>(N) + 0.01);
  }
  return ret;
}

#include <iostream>
#include <vector>
using namespace std;

int main(void) {
  int n;
  scanf("%d", &n);

  vector<double> a(n + 1), b(n + 1);
  a[0] = b[0] = 0;
  for (int i = 0; i < n; i++) {
    int tmp_a, tmp_b;
    scanf("%d %d", &tmp_a, &tmp_b);
    a[i + 1] = static_cast<double>(tmp_a);
    b[i + 1] = static_cast<double>(tmp_b);
  }
  const auto c = convolution(a, b);
  for (int i = 0; i < 2 * n; i++) {
    printf("%d\n", static_cast<int>(c[i + 1]));
  }
  return 0;
}
