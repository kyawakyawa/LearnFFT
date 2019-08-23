#include <cassert>
#include <chrono>
#include <new>
#include <random>

template <typename T>
void transpose_naive(const T* const a, const uint32_t n, const uint32_t m,
                     T* const at) {
  for (uint32_t i = 0; i < n; i++) {
    for (uint32_t j = 0; j < m; j++) {
      at[2 * (j * n + i)]     = a[2 * (i * m + j)];
      at[2 * (j * n + i) + 1] = a[2 * (i * m + j) + 1];
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

void check(const double* const a, const double* const b, const uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    assert(std::abs(a[i] - b[i]) < 1e-7);
  }
}

void test01(void) {
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::uniform_real_distribution<> dist1(0.0, 1.0);

  const size_t n1 = seed_gen() % 100 + 1;
  const size_t n2 = seed_gen() % 100 + 1;
  std::vector<double> a(n1 * n2 * 2);

  for (auto& v : a) {
    v = dist1(engine);
  }

  // for (size_t i = 0; i < n1; i++) {
  //   for (size_t j = 0; j < n2; j++) {
  //     printf("(%.3f,%.3f) ", a[2 * (i * n2 + j)], a[2 * (i * n2 + j) + 1]);
  //   }
  //   printf("\n");
  // }
  // printf("\n\n");

  std::vector<double> at(n1 * n2 * 2);
  transpose_naive(a.data(), n1, n2, at.data());

  // for (size_t i = 0; i < n2; i++) {
  //   for (size_t j = 0; j < n1; j++) {
  //     printf("(%.3f,%.3f) ", at[2 * (i * n1 + j)], at[2 * (i * n1 + j) + 1]);
  //   }
  //   printf("\n");
  // }
  // printf("\n\n");

  std::vector<double> at2(n1 * n2 * 2);
  transpose_on_cache(a.data(), n1, n2, at2.data());

  // for (size_t i = 0; i < n2; i++) {
  //   for (size_t j = 0; j < n1; j++) {
  //     printf("(%.3f,%.3f) ", at2[2 * (i * n1 + j)], at2[2 * (i * n1 + j) +
  //     1]);
  //   }
  //   printf("\n");
  // }

  check(at.data(), at2.data(), n1 * n2 * 2);
}

void test02(void) {
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::uniform_real_distribution<> dist1(0.0, 1.0);

  const size_t n1 = 2048;  // seed_gen() % 100 + 1;
  const size_t n2 = 1024;  // seed_gen() % 100 + 1;
  std::vector<double> a(n1 * n2 * 2);

  for (auto& v : a) {
    v = dist1(engine);
  }

  std::vector<double> at(n1 * n2 * 2);

  std::chrono::system_clock::time_point start, end;
  start = std::chrono::system_clock::now();
  transpose_naive(a.data(), n1, n2, at.data());
  end = std::chrono::system_clock::now();
  const double elapsed_naive =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  std::vector<double> at2(n1 * n2 * 2);

  start = std::chrono::system_clock::now();
  transpose_on_cache(a.data(), n1, n2, at2.data());
  end = std::chrono::system_clock::now();
  const double elapsed_on_cache =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  printf("naive %.10f sec   on Cache %.10f sec\n", elapsed_naive / 1000000.0,
         elapsed_on_cache / 1000000.0);
  check(at.data(), at2.data(), n1 * n2 * 2);
}
int main(void) {
  for (int i = 0; i < 100; i++) {
    test01();
  }
  for (int i = 0; i < 100; i++) {
    test02();
  }
  return 0;
}