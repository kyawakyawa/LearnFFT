#include <cmath>
#include <cstdint>
#include <vector>

#include <stdio.h>

template <typename T>
void out_with_bitrev(T* input, const uint32_t n) {
  int i = 0;
  std::vector<T> a(2 * n);
  for (size_t i = 0; i < n; i++) {
    a[i]     = input[i * 2];
    a[i + n] = input[i * 2 + 1];
  }
  T* ar = a.data();
  T* ai = a.data() + n;

  for (int j = 1; j < n - 1; j++) {
    for (int k = n >> 1; k > (i ^= k); k >>= 1)
      ;
    if (j < i) {
      const T xr = ar[j];
      const T xi = ai[j];
      ar[j]      = ar[i];
      ai[j]      = ai[i];
      ar[i]      = xr;
      ai[i]      = xi;
    }
  }
  for (size_t i = 0; i < n; i++) {
    printf("(%3.3f %3.3f) ", ar[i], ai[i]);
  }
  printf("\n");
}

template <typename T>
void cooley_turkey_fft_power2_frequency_inplace(T* a, const uint32_t n,
                                                char sign = -1) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  for (uint32_t m = n; m > 1; m >>= 1) {
    const uint32_t m2         = m / 2;
    const uint32_t num_of_fft = n / m;

    const T theta = static_cast<T>(sign) * 2.0 * M_PI / static_cast<T>(m);
    for (uint32_t i = 0; i < num_of_fft; i++) {
      for (uint32_t j = 0; j < m2; j++) {
        const uint32_t id = m * i + j;
        const T tmp_r     = a[2 * id] - a[2 * (id + m2)];
        const T tmp_i     = a[2 * id + 1] - a[2 * (id + m2) + 1];
        a[2 * id] += a[2 * (id + m2)];
        a[2 * id + 1] += a[2 * (id + m2) + 1];
        const T omega_r      = std::cos(theta * j);
        const T omega_i      = std::sin(theta * j);
        a[2 * (id + m2)]     = std::fma(tmp_r, omega_r, -tmp_i * omega_i);
        a[2 * (id + m2) + 1] = std::fma(tmp_r, omega_i, tmp_i * omega_r);
      }
    }
  }
}

template <typename T>
void cooley_turkey_fft_power2_time_inplace(T* a, const uint32_t n,
                                           char sign = -1) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  for (uint32_t m = 2; m <= n; m <<= 1) {
    const uint32_t m2         = m / 2;
    const uint32_t num_of_fft = n / m;

    const T theta = static_cast<T>(sign) * 2.0 * M_PI / static_cast<T>(m);
    for (uint32_t i = 0; i < num_of_fft; i++) {
      for (uint32_t j = 0; j < m2; j++) {
        const uint32_t id = m * i + j;
        const T omega_r   = std::cos(theta * j);
        const T omega_i   = std::sin(theta * j);
        const T tmp_r     = std::fma(a[2 * (id + m2)], omega_r,
                                 -a[2 * (id + m2) + 1] * omega_i);
        const T tmp_i =
            std::fma(a[2 * (id + m2)], omega_i, a[2 * (id + m2) + 1] * omega_r);
        a[2 * (id + m2)]     = a[2 * id] - tmp_r;
        a[2 * (id + m2) + 1] = a[2 * id + 1] - tmp_i;
        a[2 * id] += tmp_r;
        a[2 * id + 1] += tmp_i;
      }
    }
  }
}

template <typename T>
std::vector<T> convolution(const std::vector<T>& a, const std::vector<T>& b) {
  const size_t n = 2 * std::max(a.size(), b.size());
  int N          = 2;
  while (n > N) N *= 2;
  std::vector<T> ac(2 * N);
  std::vector<T> bc(2 * N);
  for (size_t i = 0; i < N; i++) {
    ac[2 * i]     = (i < a.size()) ? a[i] : static_cast<T>(0.0);
    ac[2 * i + 1] = static_cast<T>(0.0);
    bc[2 * i]     = (i < b.size()) ? b[i] : static_cast<T>(0.0);
    bc[2 * i + 1] = static_cast<T>(0.0);
  }
  cooley_turkey_fft_power2_frequency_inplace(ac.data(), N);
  cooley_turkey_fft_power2_frequency_inplace(bc.data(), N);
  for (size_t i = 0; i < N; i++) {
    const T i2    = 2 * i;
    const T tmp_r = ac[i2];
    const T tmp_i = ac[i2 + 1];
    ac[i2]        = tmp_r * bc[i2] - tmp_i * bc[i2 + 1];
    ac[i2 + 1]    = tmp_r * bc[i2 + 1] + tmp_i * bc[i2];
  }
  cooley_turkey_fft_power2_time_inplace(ac.data(), N, 1);
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
  cin >> n;

  vector<double> a(n + 1), b(n + 1);
  a[0] = b[0] = 0;
  for (int i = 0; i < n; i++) {
    cin >> a[i + 1] >> b[i + 1];
  }
  const auto c = convolution(a, b);
  for (int i = 0; i < 2 * n; i++) {
    cout << static_cast<long long>(c[i + 1]) << endl;
  }
  return 0;
}
