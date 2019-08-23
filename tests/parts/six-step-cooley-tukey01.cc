#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>

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
void transpose_on_cache(const T* const a, const uint32_t n, const uint32_t m,
                        T* const at) {
// TODO switch
#if true
  constexpr uint32_t kCacheLineSize = 64;
#else
  constexpr uint32_t kCacheLineSize =
      std::hardware_constructive_interference_size;
#endif
  constexpr uint32_t kBlock =
      std::max<uint32_t>(kCacheLineSize / sizeof(T) / 2, 1);

  for (uint32_t i = 0; i < n; i += kBlock) {
    for (uint32_t j = 0; j < m; j += kBlock) {
      const uint32_t n_ = std::min(i + kBlock, n);
      for (uint32_t i_ = i; i_ < n_; i_++) {
        const uint32_t m_ = std::min(j + kBlock, m);
        for (uint32_t j_ = j; j_ < m_; j_++) {
          at[2 * (j_ * n + i_)]     = a[2 * (i_ * m + j_)];
          at[2 * (j_ * n + i_) + 1] = a[2 * (i_ * m + j_) + 1];
        }
      }
    }
  }
}

template <typename I = uint32_t>
constexpr std::tuple<I, I> bitrev_initialization(const I n, I* const ip) {
  ip[0] = 0;
  I k   = n;
  I m   = 1;
  while (2 * m < k) {
    k = k / 2;
    for (I j = 0; j < m; j++) {
      ip[m + j] = ip[j] + k;
    }
    m = m * 2;
  }
  return std::make_tuple(m, k);
}

template <typename T, typename I = uint32_t>
constexpr void bitrv_scramble(T* const a, const I* const ip, const I n,
                              const I m, const I k) {
  if (m == k) {
    for (I i = 1; i < m; i++) {
      for (I j = 0; j < i; j++) {
        const I ji    = j + ip[i];
        const I ij    = i + ip[j];
        const T tmp_r = a[2 * ji];
        a[2 * ji]     = a[2 * ij];
        a[2 * ij]     = tmp_r;
        const T tmp_i = a[2 * ji + 1];
        a[2 * ji + 1] = a[2 * ij + 1];
        a[2 * ij + 1] = tmp_i;
      }
    }
  } else {
    for (I i = 1; i < m; i++) {
      for (I j = 0; j < i; j++) {
        const I ji     = j + ip[i];
        const I ij     = i + ip[j];
        const I tmp0_r = a[2 * ji];
        a[2 * ji]      = a[2 * ij];
        a[2 * ij]      = tmp0_r;
        const I tmp0_i = a[2 * ji + 1];
        a[2 * ji + 1]  = a[2 * ij + 1];
        a[2 * ij + 1]  = tmp0_i;

        const I tmp1_r      = a[2 * (ji + m)];
        a[2 * (ji + m)]     = a[2 * (ij + m)];
        a[2 * (ij + m)]     = tmp1_r;
        const I tmp1_i      = a[2 * (ji + m) + 1];
        a[2 * (ji + m) + 1] = a[2 * (ij + m) + 1];
        a[2 * (ij + m) + 1] = tmp1_i;
      }
    }
  }
}

template <typename T>
void twiddle(const T* const a, const uint32_t n1, const uint32_t n2,
             const T sign, T* const a_) {
  const uint32_t n = n1 * n2;

  for (uint32_t i = 0; i < n1; i++) {
    for (uint32_t j = 0; j < n2; j++) {
      const T theta     = sign * static_cast<T>(2.0) * M_PI * i * j / n;
      const T omega_r   = std::cos(theta);
      const T omega_i   = std::sin(theta);
      const uint32_t id = i * n2 + j;
      a_[2 * id]     = std::fma(a[2 * id], omega_r, -a[2 * id + 1] * omega_i);
      a_[2 * id + 1] = std::fma(a[2 * id], omega_i, a[2 * id + 1] * omega_r);
    }
  }
}

template <typename T>
void six_step_fft_power2_inplace(T* const a, const uint32_t n,
                                 const T sign = -1.0) {
  // a = {Re(a0), Im(a0), Re(a1), Im(a1), ..., Re(an-1), Im(an-1)}
  if ((n & (n - 1))) return;

  uint32_t index = 0;
  {
    uint32_t tmp = n;
    while (tmp > 1) {
      index++;
      tmp /= 2;
    }
  }
  const uint32_t index1 = index / 2;
  const uint32_t index2 = index - index1;
  const uint32_t n1     = (1 << index1);
  const uint32_t n2     = (1 << index2);

  std::vector<T> work_space(n * 2);

  transpose_on_cache(a, n1, n2, work_space.data());
  // n1r x n2c -> n2r x n1c

  for (uint32_t i = 0; i < n2; i++) {
    cooley_turkey_fft_split_radix_frequency_inplace(
        work_space.data() + 2 * n1 * i, n1, sign);
  }

  std::vector<uint32_t> ip1(n1);
  const auto [brm1, brk1] = bitrev_initialization(n1, ip1.data());
  for (uint32_t i = 0; i < n2; i++) {
    bitrv_scramble(work_space.data() + 2 * n1 * i, ip1.data(), n1, brm1, brk1);
  }

  twiddle(work_space.data(), n2, n1, sign, a);

  transpose_on_cache(a, n2, n1, work_space.data());
  // n2r x n1c -> n1r x n2c

  for (uint32_t i = 0; i < n1; i++) {
    cooley_turkey_fft_split_radix_frequency_inplace(
        work_space.data() + 2 * n2 * i, n2, sign);
  }

  std::vector<uint32_t> ip2(n2);
  const auto [brm2, brk2] = bitrev_initialization(n2, ip2.data());
  for (uint32_t i = 0; i < n1; i++) {
    bitrv_scramble(work_space.data() + 2 * n2 * i, ip2.data(), n2, brm2, brk2);
  }

  transpose_on_cache(work_space.data(), n1, n2, a);
  // n1r x n2c -> n2r x n1c
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
  six_step_fft_power2_inplace(ac.data(), N);
  six_step_fft_power2_inplace(bc.data(), N);

  for (size_t i = 0; i < N; i++) {
    const T i2    = 2 * i;
    const T tmp_r = ac[i2];
    const T tmp_i = ac[i2 + 1];
    ac[i2]        = tmp_r * bc[i2] - tmp_i * bc[i2 + 1];
    ac[i2 + 1]    = tmp_r * bc[i2 + 1] + tmp_i * bc[i2];
  }
  six_step_fft_power2_inplace(ac.data(), N, static_cast<T>(1.0));
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