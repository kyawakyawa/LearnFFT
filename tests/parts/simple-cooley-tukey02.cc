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
void cooley_turkey_fft_power4_frequency_inplace(T* a, const uint32_t n,
                                                char sign = -1) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  uint32_t m = n;
  for (; m >= 4; m >>= 2) {
    const uint32_t m4         = m / 4;
    const uint32_t num_of_fft = n / m;

    const T theta = static_cast<T>(sign) * 2.0 * M_PI / static_cast<T>(m);
    for (uint32_t i = 0; i < num_of_fft; i++) {
      for (uint32_t j = 0; j < m4; j++) {
        const uint32_t id04 = 2 * (m * i + j);
        const uint32_t id24 = id04 + 2 * (2 * m4);
        const uint32_t id14 = id04 + 2 * m4;
        const uint32_t id34 = id04 + 2 * 3 * m4;

        const T pr = a[id04] + a[id24];
        const T pi = a[id04 + 1] + a[id24 + 1];
        const T qr = a[id04] - a[id24];
        const T qi = a[id04 + 1] - a[id24 + 1];
        const T rr = a[id14] + a[id34];
        const T ri = a[id14 + 1] + a[id34 + 1];
        const T sr = -a[id14 + 1] + a[id34 + 1];
        const T si = a[id14] - a[id34];

        const T omega1r = std::cos(theta * j);
        const T omega1i = std::sin(theta * j);
        const T omega2r = std::fma(static_cast<T>(2.0) * omega1r, omega1r,
                                   -static_cast<T>(1.0));
        const T omega2i = static_cast<T>(2.0) * omega1r * omega1i;
        const T omega3r = std::fma(omega1r, omega2r, -omega1i * omega2i);
        const T omega3i = std::fma(omega1r, omega2i, omega1i * omega2r);

        a[id04]     = pr + rr;
        a[id04 + 1] = pi + ri;

        const T tmp_r1 = pr - rr;
        const T tmp_i1 = pi - ri;
        a[id14]        = std::fma(tmp_r1, omega2r, -tmp_i1 * omega2i);
        a[id14 + 1]    = std::fma(tmp_r1, omega2i, tmp_i1 * omega2r);

        const T tmp_r2 = qr - sr;
        const T tmp_i2 = qi - si;
        a[id24]        = std::fma(tmp_r2, omega1r, -tmp_i2 * omega1i);
        a[id24 + 1]    = std::fma(tmp_r2, omega1i, tmp_i2 * omega1r);

        const T tmp_r3 = qr + sr;
        const T tmp_i3 = qi + si;
        a[id34]        = std::fma(tmp_r3, omega3r, -tmp_i3 * omega3i);
        a[id34 + 1]    = std::fma(tmp_r3, omega3i, tmp_i3 * omega3r);
      }
    }
  }

  for (; m > 1; m >>= 1) {
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
  cooley_turkey_fft_power4_frequency_inplace(ac.data(), N);
  cooley_turkey_fft_power4_frequency_inplace(bc.data(), N);
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
