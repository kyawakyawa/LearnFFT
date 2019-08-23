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
constexpr void bitrv_scramble(T* const a, const T* const ip, const I n,
                              const I m, const I k) {
  if (m == k) {
    for (I i = 1; i < m; i++) {
      for (I j = 0; j < i; j++) {
        const I ji  = j + ip[i];
        const I ij  = i + ip[j];
        const T tmp = a[ji];
        a[ji]       = a[ij];
        a[ij]       = tmp;
      }
    }
  } else {
    for (I i = 1; i < m; i++) {
      for (I j = 0; j < i; j++) {
        const I ji   = j + ip[i];
        const I ij   = i + ip[j];
        const I tmp0 = a[ji];
        a[ji]        = a[ij];
        a[ij]        = tmp0;
        const I tmp1 = a[ji + m];
        a[ji + m]    = a[ij + m];
        a[ij + m]    = tmp1;
      }
    }
  }
}

int main(void) {
  constexpr uint32_t logn = 3;
  constexpr uint32_t n    = (1 << logn);

  std::vector<uint32_t> a(n), ip(n);

  for (size_t i = 0; i < n; i++) {
    a[i] = i;
  }

  const auto [m, k] = bitrev_initialization(n, ip.data());
  bitrv_scramble(a.data(), ip.data(), n, m, k);

  for (size_t i = 0; i < n; i++) {
    printf("%d\n", a[i]);
  }
  return 0;
}
