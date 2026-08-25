#include <cmath>
#include <cstddef>
#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

namespace {

// Type alias for the 3D float array
using Array3F = py::array_t<float, py::array::c_style | py::array::forcecast>;

// Compute the average and standard deviation of a 3D float array along the
// first axis, ignoring outliers beyond a specified sigma threshold. Returns a
// tuple of (avg, std, num), where avg and std are 2D arrays of shape (ny, nx)
// and num is a 2D array of shape (ny, nx) containing the number of valid pixels
// used in the computation for each (y, x) position.
py::tuple avgstd_impl(const py::array_t<float> &cube, float sigma) {
  if (cube.ndim() != 3) {
    throw std::runtime_error("cube must be 3-dimensional");
  }
  if (sigma <= 1.0f) {
    throw std::runtime_error("sigma must be greater than 1");
  }

  auto buf = cube.request();
  const std::size_t nf = static_cast<std::size_t>(buf.shape[0]);
  const std::size_t ny = static_cast<std::size_t>(buf.shape[1]);
  const std::size_t nx = static_cast<std::size_t>(buf.shape[2]);
  const float *data = static_cast<const float *>(buf.ptr);

  py::array_t<float> avg({ny, nx});
  py::array_t<float> stddev({ny, nx});
  py::array_t<std::int32_t> num({ny, nx});

  auto avg_mut = avg.mutable_unchecked<2>();
  auto std_mut = stddev.mutable_unchecked<2>();
  auto num_mut = num.mutable_unchecked<2>();

  std::vector<float> vals(nf);
  std::vector<char> ok(nf, 1);

  for (std::size_t iy = 0; iy < ny; ++iy) {
    for (std::size_t ix = 0; ix < nx; ++ix) {
      for (std::size_t iz = 0; iz < nf; ++iz) {
        vals[iz] = data[iz * ny * nx + iy * nx + ix];
        ok[iz] = 1;
      }

      std::size_t ncur = nf;
      while (true) {
        if (ncur == 0) {
          avg_mut(iy, ix) = 0.0f;
          std_mut(iy, ix) = 0.0f;
          num_mut(iy, ix) = 0;
          break;
        }

        double sum = 0.0;
        std::size_t nused = 0;
        for (std::size_t iz = 0; iz < nf; ++iz) {
          if (!ok[iz]) {
            continue;
          }
          sum += vals[iz];
          ++nused;
        }

        double tavg = nused > 0 ? sum / static_cast<double>(nused) : 0.0;
        double sumsq = 0.0;
        for (std::size_t iz = 0; iz < nf; ++iz) {
          if (!ok[iz]) {
            continue;
          }
          const double diff = vals[iz] - tavg;
          sumsq += diff * diff;
        }

        double tstd =
            nused > 1 ? std::sqrt(sumsq / static_cast<double>(nused - 1)) : 0.0;
        const double thresh = sigma * tstd;

        std::vector<char> new_ok(ok.begin(), ok.end());
        std::size_t nnew = 0;
        for (std::size_t iz = 0; iz < nf; ++iz) {
          const bool keep = ok[iz] && std::fabs(vals[iz] - tavg) <= thresh;
          new_ok[iz] = keep ? 1 : 0;
          if (keep) {
            ++nnew;
          }
        }

        if (nnew == ncur) {
          double final_sum = 0.0;
          for (std::size_t iz = 0; iz < nf; ++iz) {
            if (!new_ok[iz]) {
              continue;
            }
            final_sum += vals[iz];
          }
          const double final_avg =
              nnew > 0 ? final_sum / static_cast<double>(nnew) : 0.0;
          double final_sumsq = 0.0;
          for (std::size_t iz = 0; iz < nf; ++iz) {
            if (!new_ok[iz]) {
              continue;
            }
            const double diff = vals[iz] - final_avg;
            final_sumsq += diff * diff;
          }
          const double final_std =
              nnew > 1 ? std::sqrt(final_sumsq / static_cast<double>(nnew - 1))
                       : 0.0;
          avg_mut(iy, ix) = static_cast<float>(final_avg);
          std_mut(iy, ix) = static_cast<float>(final_std);
          num_mut(iy, ix) = static_cast<std::int32_t>(nnew);
          break;
        }

        ok.swap(new_ok);
        ncur = nnew;
      }
    }
  }

  return py::make_tuple(avg, stddev, num);
}

} // namespace

PYBIND11_MODULE(_support_cpp, m) {
  m.def("avgstd", &avgstd_impl, py::arg("cube"), py::arg("sigma"));
}
