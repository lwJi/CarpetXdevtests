#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "loop.hxx"

#define GFINDEXTYPE 1, 1, 1

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t, const T x, const T y,
                        const T z, T &u, T &rho) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    rho = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
  } else {
    u = (f(t - r) - f(t + r)) / r;
    rho = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
  }
}

extern "C" void TestRestrict_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestRestrict_Init;

  CCTK_INFO("Initializing grid function");
  grid.loop_all<GFINDEXTYPE>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    CCTK_REAL u, rho;
    gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u, rho);
    iteration(p.I) = u;
  });
}

extern "C" void TestRestrict_Update(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestRestrict_Update;

  CCTK_VINFO("Updating grid function at iteration %d level %d time %g",
             cctk_iteration, cctk_level, cctk_time);
  grid.loop_int<GFINDEXTYPE>(
      grid.nghostzones, [=] CCTK_HOST(const Loop::PointDesc &p) {
        CCTK_REAL u, rho;
        gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u, rho);
        iteration(p.I) = u + cctk_level;
      });
}
