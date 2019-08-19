#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

template <typename T>
static std::vector<std::vector<std::complex<T>>> omega;
template <typename T>
static std::vector<std::vector<std::complex<T>>> omega_inv;
template <typename T>
void make_omega_table(const uint_fast8_t logn) {
  size_t i = omega<T>.size();
  if (i >= logn) return;

  for (; i <= logn; i++) {
    omega<T>.emplace_back(std::vector<std::complex<T>>());
    omega_inv<T>.emplace_back(std::vector<std::complex<T>>());
    const uint32_t n = (1 << i);
    const T invn     = static_cast<T>(1.0) / n;
    for (size_t j = 0; j < n / 2; j++) {
      const T tmp = static_cast<T>(-2.0) * M_PI * j * invn;
      omega<T>.back().emplace_back(std::cos(tmp), std::sin(tmp));
      omega_inv<T>.back().emplace_back(std::cos(-tmp), std::sin(-tmp));
    }
  }
}

static std::vector<std::vector<uint32_t>> bitrev;
void make_bitrev_table(const uint_fast8_t logn) {
  size_t i = bitrev.size();
  if (i >= logn) return;

  for (; i <= logn; i++) {
    bitrev.emplace_back(std::vector<uint32_t>());
    const uint32_t n = (1 << i);
    for (uint32_t j = 0; j < n; j++) bitrev.back().emplace_back(j);

    std::transform(bitrev.back().begin(), bitrev.back().end(),
                   bitrev.back().begin(), [&i](uint32_t t) {
                     uint32_t ret = 0;
                     for (uint32_t j = 0; j < i; j++) {
                       ret <<= 1;
                       ret |= (t & 1);
                       t >>= 1;
                     }
                     return ret;
                   });
  }
}

template <typename T>
static void cooley_tukey_cfft_frequency_inplace(std::complex<T> *a,
                                                const uint_fast8_t logn) {
  const uint32_t n = (1 << logn);
  for (uint32_t i = 0; i < logn; i++) {
    const uint32_t offset = (n >> (i + 1));
    const uint32_t m      = (1 << i);
    for (uint32_t j = 0; j < m; j++) {
      for (uint32_t k = 0; k < offset; k++) {
        const uint32_t id = j * 2 * offset + k;
        const auto &tmp_a = a[id];
        const auto &tmp_b = a[id + offset];
        const auto tmp    = tmp_a - tmp_b;
        a[id]             = tmp_a + tmp_b;
        a[id + offset]    = tmp * omega<T>[logn - i][k];
      }
    }
  }
}

template <typename T>
static void cooley_tukey_cfft_time_inplace(std::complex<T> *a,
                                           const uint_fast8_t logn) {
  const uint32_t n = (1 << logn);
  for (uint32_t i = 0; i < logn; i++) {
    const uint32_t offset = (1 << i);
    const uint32_t m      = (n >> (i + 1));
    for (uint32_t j = 0; j < m; j++) {
      for (uint32_t k = 0; k < offset; k++) {
        const uint32_t id = j * 2 * offset + k;
        a[id + offset] *= omega<T>[i + 1][k];
        const auto &tmp_a = a[id];
        const auto &tmp_b = a[id + offset];
        const auto tmp    = tmp_a - tmp_b;
        a[id]             = tmp_a + tmp_b;
        a[id + offset]    = tmp;
      }
    }
  }
}

template <typename T>
void cfft_inplace(std::vector<std::complex<T>> *a) {
  if (a->empty()) return;
  const size_t n = a->size();

  uint_fast8_t logn = 0;
  while (!((1 << logn) & n)) logn++;

  make_omega_table<T>(logn);
  make_bitrev_table(logn);

  cooley_tukey_cfft_frequency_inplace(a->data(), logn);
  for (size_t i = 0; i < n; i++) {
    if (i < bitrev[logn][i])
      swap(a->operator[](i), a->operator[](bitrev[logn][i]));
  }
}

template <typename T>
static void cooley_tukey_cinvfft_frequency_inplace(std::complex<T> *a,
                                                   const uint_fast8_t logn) {
  const uint32_t n = (1 << logn);
  for (uint32_t i = 0; i < logn; i++) {
    const uint32_t offset = (n >> (i + 1));
    const uint32_t m      = (1 << i);
    for (uint32_t j = 0; j < m; j++) {
      for (uint32_t k = 0; k < offset; k++) {
        const uint32_t id = j * 2 * offset + k;
        const auto &tmp_a = a[id];
        const auto &tmp_b = a[id + offset];
        const auto tmp    = tmp_a - tmp_b;
        a[id]             = tmp_a + tmp_b;
        a[id + offset]    = tmp * omega_inv<T>[logn - i][k];
      }
    }
  }

  for (uint32_t i = 0; i < n; i++) a[i] /= n;
}

template <typename T>
static void cooley_tukey_cinvfft_time_inplace(std::complex<T> *a,
                                              const uint8_t logn) {
  const uint32_t n = (1 << logn);
  for (uint32_t i = 0; i < logn; i++) {
    const uint32_t offset = (1 << i);
    const uint32_t m      = (n >> (i + 1));
    for (uint32_t j = 0; j < m; j++) {
      for (uint32_t k = 0; k < offset; k++) {
        const uint32_t id = j * 2 * offset + k;
        a[id + offset] *= omega_inv<T>[i + 1][k];
        const auto &tmp_a = a[id];
        const auto &tmp_b = a[id + offset];
        const auto tmp    = tmp_a - tmp_b;
        a[id]             = tmp_a + tmp_b;
        a[id + offset]    = tmp;
      }
    }
  }

  for (uint32_t i = 0; i < n; i++) a[i] /= n;
}

template <typename T>
void cinvfft_inplace(std::vector<std::complex<T>> *a) {
  if (a->empty()) return;
  const size_t n = a->size();

  uint_fast8_t logn = 0;
  while (!((1 << logn) & n)) logn++;

  make_omega_table<T>(logn);
  make_bitrev_table(logn);

  for (size_t i = 0; i < n; i++) {
    if (i < bitrev[logn][i])
      swap(a->operator[](i), a->operator[](bitrev[logn][i]));
  }

  cooley_tukey_cinvfft_time_inplace(a->data(), logn);
}

#include <iostream>

int convolution(long long a[], long long b[], long long c[], int n) {
  n *= 2;
  int N = 2, index = 1;
  ;
  while (n > N) {
    N *= 2;
    index++;
  }
  std::vector<std::complex<double>> ac(N);
  std::vector<std::complex<double>> bc(N);

  for (int i = 0; i < N; i++) {
    ac[i] = std::complex<double>((i < n / 2) ? a[i] : 0.0, 0.0);
    bc[i] = std::complex<double>((i < n / 2) ? b[i] : 0.0, 0.0);
  }
  cfft_inplace<double>(&ac);
  cfft_inplace<double>(&bc);
  for (int i = 0; i < N; i++) {
    ac.at(i) *= bc.at(i);
  }
  cinvfft_inplace<double>(&ac);
  for (int i = 0; i < n; i++) c[i] = ac[i].real() + 0.01;
  return N;
}

int convolution2(long long a[], long long b[], long long c[], int n) {
  n *= 2;
  int N = 2, index = 1;
  ;
  while (n > N) {
    N *= 2;
    index++;
  }
  std::vector<std::complex<double>> ac(N);

  for (int i = 0; i < N; i++) {
    ac[i] = std::complex<double>((i < n / 2) ? a[i] : 0.0,
                                 (i < n / 2) ? b[i] : 0.0);
  }
  cfft_inplace<double>(&ac);
  ac[0]     = ac[0].real() * ac[0].imag();
  ac[N / 2] = ac[N / 2].real() * ac[N / 2].imag();
  for (int i = 1; i < N / 2; i++) {
    const auto &&tmp = ac[i] * ac[i] - std::conj(ac[N - i] * ac[N - i]);
    ac[i] = std::complex<double>(0.25 * tmp.imag(), -0.25 * tmp.real());
  }
  for (int i = N / 2; i < N; i++) {
    ac[i] = std::conj(ac[N - i]);
  }

  {
    const auto tmp_r = (ac[0].real() + ac[N / 2].real()) / 2.0;
    const auto tmp_i = (ac[0].real() - ac[N / 2].real()) / 2.0;

    ac[0] = std::complex<double>(tmp_r, tmp_i);
  }
  for (int i = 1; i < N / 4 /*floor*/; i++) {
    const auto tmp =
        (std::complex<double>(0.0, 1.0) * std::conj(omega<double>[index][i]) -
         1.0) *
        0.5;
    ac[N / 2 - i]      = std::conj(ac[N / 2 - i]);
    const auto &&tmp_c = tmp * (ac[i] - ac[N / 2 - i]);

    ac[i] += tmp_c;
    ac[N / 2 - i] -= tmp_c;
    ac[N / 2 - i] = std::conj(ac[N / 2 - i]);
  }
  if (N % 4 == 0) {
    ac[N / 4] = std::conj(ac[N / 4]);
  }
  ac.resize(N / 2);
  cinvfft_inplace<double>(&ac);
  for (int i = 0; i < n / 2; i++) {
    c[i * 2 + 0] = ac[i].real() + 0.01;
    c[i * 2 + 1] = ac[i].imag() + 0.01;
  }
  return N;
}

long long n;
long long a[100001], b[100001], c[200002];

int main() {
  std::cin >> n;

  for (int i = 0; i < n; i++)
    std::cin >> a[i + 1] >> b[i + 1];  // a[0] = b[0] = 0とする

  convolution2(a, b, c, n + 1);  // 0〜nのn+1個
  for (int i = 0; i < 2 * n; i++) std::cout << c[i + 1] << std::endl;
}
