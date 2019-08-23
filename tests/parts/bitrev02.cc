#include <tuple>
#include <vector>

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

int main(void) {
  constexpr uint32_t logn = 3;
  constexpr uint32_t n    = (1 << logn);

  std::vector<uint32_t> a(n * 2), ip(n);

  for (size_t i = 0; i < n; i++) {
    a[2 * i]     = i;
    a[2 * i + 1] = i;
  }

  const auto [m, k] = bitrev_initialization(n, ip.data());
  bitrv_scramble(a.data(), ip.data(), n, m, k);

  for (size_t i = 0; i < 2 * n; i++) {
    printf("%d\n", a[i]);
  }
  return 0;
}
