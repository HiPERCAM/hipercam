#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

// Portable wrapper for compiler-specific restrict support
#if defined(__clang__) || defined(__GNUC__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace {

// Helper function to generate sub-pixel offsets for binning.
// This creates a vector of offsets for each bin and sub-bin, centered around
// zero, to be used in the sub-pixellation loops in the profile evaluation and
// derivatives.
std::vector<double> make_subpixel_offsets(int bin, int ndiv) {
  std::vector<double> offsets;
  if (ndiv <= 0) {
    return offsets;
  }
  offsets.reserve(static_cast<std::size_t>(bin * ndiv));
  double inv_ndiv = 1.0 / static_cast<double>(ndiv);
  double soff = (ndiv - 1.0) / (2.0 * ndiv);

  for (int ibin = 0; ibin < bin; ++ibin) {
    double base_off = ibin - (bin - 1) / 2.0 - soff;
    for (int isub = 0; isub < ndiv; ++isub) {
      offsets.push_back(base_off + isub * inv_ndiv);
    }
  }

  return offsets;
}

// Helper function to calculate the Moffat alpha parameter.
double calc_moffat_alpha(double fwhm, double beta) {
  double tbeta = std::max(0.01, beta);
  return 4.0 * (std::pow(2.0, 1.0 / tbeta) - 1.0) / (fwhm * fwhm);
}

// Evaluate one Moffat model value at a single pixel coordinate.
// This is used by the selected-pixel residual path to avoid building
// full 2D model arrays when only masked pixels are needed.
double moffat_value_at(double x_val, double y_val, double height, double xcen,
                       double ycen, double alpha, double tbeta,
                       const std::vector<double> &x_offsets,
                       const std::vector<double> &y_offsets) {
  if (!x_offsets.empty() && !y_offsets.empty()) {
    double prof = 0.0;
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

    for (double yoff : y_offsets) {
      for (double xoff : x_offsets) {
        double dx = x_val + xoff - xcen;
        double dy = y_val + yoff - ycen;
        double rsq = dx * dx + dy * dy;
        prof += std::pow(1.0 + alpha * rsq, -tbeta);
      }
    }
    return height * inv_nadd * prof;
  }

  double dx = x_val - xcen;
  double dy = y_val - ycen;
  double rsq = dx * dx + dy * dy;
  return height * std::pow(1.0 + alpha * rsq, -tbeta);
}

// Evaluate one Gaussian model value at a single pixel coordinate.
// This is used by the selected-pixel residual path to avoid building
// full 2D model arrays when only masked pixels are needed.
double gaussian_value_at(double x_val, double y_val, double height, double xcen,
                         double ycen, double alpha,
                         const std::vector<double> &x_offsets,
                         const std::vector<double> &y_offsets) {
  if (!x_offsets.empty() && !y_offsets.empty()) {
    double prof = 0.0;
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

    for (double yoff : y_offsets) {
      for (double xoff : x_offsets) {
        double dx = x_val + xoff - xcen;
        double dy = y_val + yoff - ycen;
        double rsq = dx * dx + dy * dy;
        prof += std::exp(-alpha * rsq);
      }
    }
    return height * inv_nadd * prof;
  }

  double dx = x_val - xcen;
  double dy = y_val - ycen;
  double rsq = dx * dx + dy * dy;
  return height * std::exp(-alpha * rsq);
}

// Calculate Moffat derivatives at one pixel coordinate.
// The outputs follow the same normalization as dmoffat so the fit path
// stays numerically equivalent.
void moffat_derivs_at(double x_val, double y_val, double height, double xcen,
                      double ycen, double alpha, double tbeta,
                      double dfwhm_coeff, double dbeta_coeff,
                      const std::vector<double> &x_offsets,
                      const std::vector<double> &y_offsets, bool comp_dfwhm,
                      bool comp_dbeta, double &dheight, double &dxcen,
                      double &dycen, double &dfwhm, double &dbeta) {
  double two_alpha_tbeta = 2.0 * alpha * tbeta;
  dheight = 0.0;
  dxcen = 0.0;
  dycen = 0.0;
  dfwhm = 0.0;
  dbeta = 0.0;

  if (!x_offsets.empty() && !y_offsets.empty()) {
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

    // Instead of checking comp_dfwhm and comp_dbeta inside the innermost loop,
    // we split into four separate loops here to save time on the checks.
    if (comp_dfwhm && comp_dbeta) {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double denom = 1.0 + alpha * rsq;
          // Compute denom^-tbeta once, then derive save1 from it to avoid
          // paying for a second pow call in the loop.
          double dh = std::pow(denom, -tbeta);
          double save1 = height * dh / denom;
          double save2 = save1 * rsq;
          dheight += dh;
          dxcen += two_alpha_tbeta * dx * save1;
          dycen += two_alpha_tbeta * dy * save1;
          dfwhm += dfwhm_coeff * save2;
          double log_denom = std::log(denom);
          dbeta += (-log_denom * height * dh + dbeta_coeff * save2);
        }
      }
    } else if (comp_dfwhm) {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double denom = 1.0 + alpha * rsq;
          double dh = std::pow(denom, -tbeta);
          double save1 = height * dh / denom;
          double save2 = save1 * rsq;
          dheight += dh;
          dxcen += two_alpha_tbeta * dx * save1;
          dycen += two_alpha_tbeta * dy * save1;
          dfwhm += dfwhm_coeff * save2;
        }
      }
    } else if (comp_dbeta) {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double denom = 1.0 + alpha * rsq;
          double dh = std::pow(denom, -tbeta);
          double save1 = height * dh / denom;
          double save2 = save1 * rsq;
          dheight += dh;
          dxcen += two_alpha_tbeta * dx * save1;
          dycen += two_alpha_tbeta * dy * save1;
          double log_denom = std::log(denom);
          dbeta += (-log_denom * height * dh + dbeta_coeff * save2);
        }
      }
    } else {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double denom = 1.0 + alpha * rsq;
          double dh = std::pow(denom, -tbeta);
          double save1 = height * dh / denom;
          dheight += dh;
          dxcen += two_alpha_tbeta * dx * save1;
          dycen += two_alpha_tbeta * dy * save1;
        }
      }
    }

    dheight *= inv_nadd;
    dxcen *= inv_nadd;
    dycen *= inv_nadd;
    if (comp_dfwhm) {
      dfwhm *= inv_nadd;
    }
    if (comp_dbeta) {
      dbeta *= inv_nadd;
    }
    return;
  }

  double dx = x_val - xcen;
  double dy = y_val - ycen;
  double rsq = dx * dx + dy * dy;
  double denom = 1.0 + alpha * rsq;
  dheight = std::pow(denom, -tbeta);
  double save1 = height * dheight / denom;
  double save2 = save1 * rsq;
  dxcen = two_alpha_tbeta * dx * save1;
  dycen = two_alpha_tbeta * dy * save1;
  if (comp_dfwhm) {
    dfwhm = dfwhm_coeff * save2;
  }
  if (comp_dbeta) {
    double log_denom = std::log(denom);
    dbeta = (-log_denom * height * dheight + dbeta_coeff * save2);
  }
}

// Calculate Gaussian derivatives at one pixel coordinate.
// The outputs follow the same normalization as dgaussian so the fit path
// stays numerically equivalent.
void gaussian_derivs_at(double x_val, double y_val, double height, double xcen,
                        double ycen, double alpha, double two_alpha_height,
                        double dfwhm_coeff,
                        const std::vector<double> &x_offsets,
                        const std::vector<double> &y_offsets, bool comp_dfwhm,
                        double &dheight, double &dxcen, double &dycen,
                        double &dfwhm) {
  dheight = 0.0;
  dxcen = 0.0;
  dycen = 0.0;
  dfwhm = 0.0;

  if (!x_offsets.empty() && !y_offsets.empty()) {
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

    // Same comp_dfwhm check optimization as in moffat_derivs_at,
    // except here we don't have comp_dbeta to worry about.
    if (comp_dfwhm) {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double dh = std::exp(-alpha * rsq);
          dheight += dh;
          dxcen += two_alpha_height * dh * dx;
          dycen += two_alpha_height * dh * dy;
          dfwhm += dfwhm_coeff * dh * rsq;
        }
      }
    } else {
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          double dh = std::exp(-alpha * rsq);
          dheight += dh;
          dxcen += two_alpha_height * dh * dx;
          dycen += two_alpha_height * dh * dy;
        }
      }
    }

    dheight *= inv_nadd;
    dxcen *= inv_nadd;
    dycen *= inv_nadd;
    if (comp_dfwhm) {
      dfwhm *= inv_nadd;
    }
    return;
  }

  double dx = x_val - xcen;
  double dy = y_val - ycen;
  double rsq = dx * dx + dy * dy;

  dheight = std::exp(-alpha * rsq);
  dxcen = two_alpha_height * dheight * dx;
  dycen = two_alpha_height * dheight * dy;
  if (comp_dfwhm) {
    dfwhm = dfwhm_coeff * dheight * rsq;
  }
}

// Shared shape validation for profile, residual, and Jacobian entry points.
// py::buffer_info is move-only, so callers pass pointers rather than values.
void check_matching_shapes(
    std::initializer_list<const py::buffer_info *> infos) {
  const py::buffer_info &first = **infos.begin();
  for (const py::buffer_info *info : infos) {
    if (info->ndim != first.ndim || info->shape != first.shape) {
      throw std::runtime_error("Input arrays must have the same shape");
    }
  }
}

// Companion check for the 1D ok_indices array passed with shaped inputs.
void check_ok_indices_1d(const py::buffer_info &ok_info) {
  if (ok_info.ndim != 1) {
    throw std::runtime_error(
        "ok_indices must be a 1D array");
  }
}

// Bounds-check every selected pixel index against the flattened pixel count.
void check_ok_indices_range(const std::int64_t *ok_ptr, std::size_t n_ok,
                            std::size_t n_pixels) {
  for (std::size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];
    if (idx < 0 || static_cast<std::size_t>(idx) >= n_pixels) {
      throw std::runtime_error("ok_indices contains out-of-range values");
    }
  }
}

} // namespace

py::array_t<double> moffat_resid(py::array_t<double> x, py::array_t<double> y,
                                 py::array_t<double> data,
                                 py::array_t<double> sigma,
                                 py::array_t<std::int64_t> ok_indices,
                                 double sky, double height, double xcen,
                                 double ycen, double fwhm, double beta,
                                 int xbin, int ybin, int ndiv) {

  // Compute residuals only at valid indices provided by Python. This
  // bypasses full-frame residual assembly and boolean masking in Python.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info data_info = data.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  check_matching_shapes({&x_info, &y_info, &data_info, &sigma_info});
  check_ok_indices_1d(ok_info);

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *data_ptr = static_cast<double *>(data_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);
  std::size_t n_ok = ok_info.shape[0];

  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(static_cast<py::ssize_t>(n_ok));
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  check_ok_indices_range(ok_ptr, n_ok, n_pixels);

  // Release the GIL now the Python operations are complete.
  // pybind11 buffer operations (request, array creation) require the GIL,
  // so we can only release it after extracting all pointers and dimensions.
  py::gil_scoped_release release;

  for (std::size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double model =
        sky + moffat_value_at(x_ptr[idx], y_ptr[idx], height, xcen, ycen, alpha,
                              tbeta, x_offsets, y_offsets);
    result_ptr[i] = (data_ptr[idx] - model) / sigma_ptr[idx];
  }

  return result;
}

py::array_t<double> dmoffat_jac(py::array_t<double> x, py::array_t<double> y,
                                py::array_t<double> sigma,
                                py::array_t<std::int64_t> ok_indices,
                                double sky, double height, double xcen,
                                double ycen, double fwhm, double beta, int xbin,
                                int ybin, int ndiv, bool comp_dfwhm,
                                bool comp_dbeta, const std::vector<int> &inds) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Build Jacobian rows directly for selected pixels and selected parameter
  // columns (`inds`), matching the same derivative ordering used by Python.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  check_matching_shapes({&x_info, &y_info, &sigma_info});
  check_ok_indices_1d(ok_info);

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);
  std::size_t n_ok = ok_info.shape[0];
  std::size_t n_par = inds.size();

  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);
  double two_alpha_tbeta = 2.0 * alpha * tbeta;
  double dfwhm_coeff = two_alpha_tbeta / fwhm;
  double dbeta_coeff =
      4.0 * std::log(2.0) * std::pow(2.0, 1.0 / tbeta) / tbeta / (fwhm * fwhm);
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(
      {static_cast<py::ssize_t>(n_ok), static_cast<py::ssize_t>(n_par)});
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  check_ok_indices_range(ok_ptr, n_ok, n_pixels);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (std::size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double dheight, dxcen, dycen, dfwhm, dbeta;
    moffat_derivs_at(x_ptr[idx], y_ptr[idx], height, xcen, ycen, alpha, tbeta,
                     dfwhm_coeff, dbeta_coeff, x_offsets, y_offsets, comp_dfwhm,
                     comp_dbeta, dheight, dxcen, dycen, dfwhm, dbeta);

    double d0 = 1.0;
    double d1 = dheight;
    double d2 = dxcen;
    double d3 = dycen;
    double d4;
    double d5;

    // Keep compatibility with the legacy dmoffat API that duplicates
    // placeholders when derivatives are not requested.
    if (comp_dfwhm && comp_dbeta) {
      d4 = dfwhm;
      d5 = dbeta;
    } else if (comp_dfwhm) {
      d4 = dfwhm;
      d5 = dfwhm;
    } else if (comp_dbeta) {
      d4 = dbeta;
      d5 = dbeta;
    } else {
      d4 = dycen;
      d5 = dycen;
    }

    double derivs[6] = {d0, d1, d2, d3, d4, d5};
    double inv_sigma = -1.0 / sigma_ptr[idx];

    for (std::size_t j = 0; j < n_par; ++j) {
      int ind = inds[j];
      if (ind < 0 || ind > 5) {
        throw std::runtime_error("inds contains out-of-range derivative index");
      }
      result_ptr[i * n_par + j] = derivs[ind] * inv_sigma;
    }
  }

  return result;
}

py::array_t<double> gaussian_resid(py::array_t<double> x, py::array_t<double> y,
                                   py::array_t<double> data,
                                   py::array_t<double> sigma,
                                   py::array_t<std::int64_t> ok_indices,
                                   double sky, double height, double xcen,
                                   double ycen, double fwhm, int xbin, int ybin,
                                   int ndiv) {

  // Gaussian equivalent of moffat_resid: selected-pixel residuals only.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info data_info = data.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  check_matching_shapes({&x_info, &y_info, &data_info, &sigma_info});
  check_ok_indices_1d(ok_info);

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *data_ptr = static_cast<double *>(data_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);
  std::size_t n_ok = ok_info.shape[0];

  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(static_cast<py::ssize_t>(n_ok));
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  check_ok_indices_range(ok_ptr, n_ok, n_pixels);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (std::size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double model = sky + gaussian_value_at(x_ptr[idx], y_ptr[idx], height, xcen,
                                           ycen, alpha, x_offsets, y_offsets);
    result_ptr[i] = (data_ptr[idx] - model) / sigma_ptr[idx];
  }

  return result;
}

py::array_t<double> dgaussian_jac(py::array_t<double> x, py::array_t<double> y,
                                  py::array_t<double> sigma,
                                  py::array_t<std::int64_t> ok_indices,
                                  double sky, double height, double xcen,
                                  double ycen, double fwhm, int xbin, int ybin,
                                  int ndiv, bool comp_dfwhm,
                                  const std::vector<int> &inds) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Selected-pixel Jacobian assembly for Gaussian fits.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  check_matching_shapes({&x_info, &y_info, &sigma_info});
  check_ok_indices_1d(ok_info);

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);
  std::size_t n_ok = ok_info.shape[0];
  std::size_t n_par = inds.size();

  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  double two_alpha_height = 2.0 * alpha * height;
  double dfwhm_coeff = two_alpha_height / fwhm;
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(
      {static_cast<py::ssize_t>(n_ok), static_cast<py::ssize_t>(n_par)});
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  check_ok_indices_range(ok_ptr, n_ok, n_pixels);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (std::size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double dheight, dxcen, dycen, dfwhm;
    gaussian_derivs_at(x_ptr[idx], y_ptr[idx], height, xcen, ycen, alpha,
                       two_alpha_height, dfwhm_coeff, x_offsets, y_offsets,
                       comp_dfwhm, dheight, dxcen, dycen, dfwhm);

    double d0 = 1.0;
    double d1 = dheight;
    double d2 = dxcen;
    double d3 = dycen;
    double d4 = comp_dfwhm ? dfwhm : dycen;

    double derivs[5] = {d0, d1, d2, d3, d4};
    double inv_sigma = -1.0 / sigma_ptr[idx];

    for (std::size_t j = 0; j < n_par; ++j) {
      int ind = inds[j];
      if (ind < 0 || ind > 4) {
        throw std::runtime_error("inds contains out-of-range derivative index");
      }
      result_ptr[i * n_par + j] = derivs[ind] * inv_sigma;
    }
  }

  return result;
}

// C++ implementation of the Moffat profile function
py::array_t<double> moffat(py::array_t<double> x, py::array_t<double> y,
                           double sky, double height, double xcen, double ycen,
                           double fwhm, double beta, int xbin, int ybin,
                           int ndiv) {

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  check_matching_shapes({&x_info, &y_info});

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output array with same shape as input
  py::array_t<double> result = py::array_t<double>(x_info.shape);
  py::buffer_info result_info = result.request();
  double *RESTRICT result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Moffat profile parameters
  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  if (ndiv > 0) {
    // With sub-pixellation
    const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
    const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

// OpenMP parallelization of for loop with SIMD vectorization.
#ifdef _OPENMP
#pragma omp parallel for simd
#endif

    // Loop over all pixels
    for (std::size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double prof = 0.0;

      // Reuse precomputed sub-pixel offsets across all pixels.
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          prof += std::pow(1.0 + alpha * rsq, -tbeta);
        }
      }

      result_ptr[pixel_idx] = sky + height * inv_nadd * prof;
    }
  } else {
    // Fast calculation at pixel centers
    for (std::size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;
      result_ptr[i] = sky + height * std::pow(1.0 + alpha * rsq, -tbeta);
    }
  }

  return result;
}

// C++ implementation of the Moffat profile derivatives
std::vector<py::array_t<double>>
dmoffat(py::array_t<double> x, py::array_t<double> y, double sky, double height,
        double xcen, double ycen, double fwhm, double beta, int xbin, int ybin,
        int ndiv, bool comp_dfwhm, bool comp_dbeta) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  check_matching_shapes({&x_info, &y_info});

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output arrays with same shape as input
  std::vector<py::array_t<double>> result;
  result.reserve(6);

  // Always need dsky, dheight, dxcen, dycen
  py::array_t<double> dsky = py::array_t<double>(x_info.shape);
  py::array_t<double> dheight = py::array_t<double>(x_info.shape);
  py::array_t<double> dxcen = py::array_t<double>(x_info.shape);
  py::array_t<double> dycen = py::array_t<double>(x_info.shape);

  // Optionally need dfwhm and dbeta
  py::array_t<double> dfwhm;
  py::array_t<double> dbeta;

  if (comp_dfwhm) {
    dfwhm = py::array_t<double>(x_info.shape);
  }

  if (comp_dbeta) {
    dbeta = py::array_t<double>(x_info.shape);
  }

  // Get data pointers
  py::buffer_info dsky_info = dsky.request();
  py::buffer_info dheight_info = dheight.request();
  py::buffer_info dxcen_info = dxcen.request();
  py::buffer_info dycen_info = dycen.request();

  double *RESTRICT dsky_ptr = static_cast<double *>(dsky_info.ptr);
  double *RESTRICT dheight_ptr = static_cast<double *>(dheight_info.ptr);
  double *RESTRICT dxcen_ptr = static_cast<double *>(dxcen_info.ptr);
  double *RESTRICT dycen_ptr = static_cast<double *>(dycen_info.ptr);

  py::buffer_info dfwhm_info;
  py::buffer_info dbeta_info;
  double *RESTRICT dfwhm_ptr = nullptr;
  double *RESTRICT dbeta_ptr = nullptr;

  if (comp_dfwhm) {
    dfwhm_info = dfwhm.request();
    dfwhm_ptr = static_cast<double *>(dfwhm_info.ptr);
  }

  if (comp_dbeta) {
    dbeta_info = dbeta.request();
    dbeta_ptr = static_cast<double *>(dbeta_info.ptr);
  }

  // Calculate Moffat profile parameters
  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);
  double two_alpha_tbeta = 2.0 * alpha * tbeta;
  double dfwhm_coeff = two_alpha_tbeta / fwhm;
  double dbeta_coeff =
      4.0 * std::log(2.0) * std::pow(2.0, 1.0 / tbeta) / tbeta / (fwhm * fwhm);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  // Initialize dsky to ones (derivative of sky is always 1)
  std::fill_n(dsky_ptr, n_pixels, 1.0);

  if (ndiv > 0) {
    // With sub-pixellation
    const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
    const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

// OpenMP parallelization of for loop with SIMD vectorization.
#ifdef _OPENMP
#pragma omp parallel for simd
#endif

    // Loop over all pixels
    for (std::size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      // Use _sum suffix to avoid shadowing the outer py::array_t declarations.
      double dheight_sum = 0.0;
      double dxcen_sum = 0.0;
      double dycen_sum = 0.0;
      double dfwhm_sum = 0.0;
      double dbeta_sum = 0.0;

      // Reuse precomputed sub-pixel offsets across all pixels.
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;

          double denom = 1.0 + alpha * rsq;
          // Derivatives
          // Keep the ordering here aligned with moffat_derivs_at: compute
          // denom^-tbeta once, then derive save1 from it to avoid a second pow.
          double dh = std::pow(denom, -tbeta);
          double save1 = height * dh / denom;
          double save2 = save1 * rsq;
          dheight_sum += dh;
          dxcen_sum += two_alpha_tbeta * dx * save1;
          dycen_sum += two_alpha_tbeta * dy * save1;

          if (comp_dfwhm) {
            dfwhm_sum += dfwhm_coeff * save2;
          }

          if (comp_dbeta) {
            double log_denom = std::log(denom);
            dbeta_sum += (-log_denom * height * dh + dbeta_coeff * save2);
          }
        }
      }

      dheight_ptr[pixel_idx] = dheight_sum * inv_nadd;
      dxcen_ptr[pixel_idx] = dxcen_sum * inv_nadd;
      dycen_ptr[pixel_idx] = dycen_sum * inv_nadd;

      if (comp_dfwhm) {
        dfwhm_ptr[pixel_idx] = dfwhm_sum * inv_nadd;
      }

      if (comp_dbeta) {
        dbeta_ptr[pixel_idx] = dbeta_sum * inv_nadd;
      }
    }
  } else {
    // Fast calculation at pixel centers
    for (std::size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;

      double denom = 1.0 + alpha * rsq;
      dheight_ptr[i] = std::pow(denom, -tbeta);
      double save1 = height * dheight_ptr[i] / denom;
      double save2 = save1 * rsq;

      // Derivatives
      dxcen_ptr[i] = two_alpha_tbeta * dx * save1;
      dycen_ptr[i] = two_alpha_tbeta * dy * save1;

      if (comp_dfwhm) {
        dfwhm_ptr[i] = dfwhm_coeff * save2;
      }

      if (comp_dbeta) {
        double log_denom = std::log(denom);
        dbeta_ptr[i] =
            (-log_denom * height * dheight_ptr[i] + dbeta_coeff * save2);
      }
    }
  }

  // Add arrays to result vector
  result.push_back(dsky);
  result.push_back(dheight);
  result.push_back(dxcen);
  result.push_back(dycen);

  if (comp_dfwhm && comp_dbeta) {
    result.push_back(dfwhm);
    result.push_back(dbeta);
  } else if (comp_dfwhm) {
    result.push_back(dfwhm);
    result.push_back(dfwhm); // duplicate to match Python API
  } else if (comp_dbeta) {
    result.push_back(dbeta);
    result.push_back(dbeta); // duplicate to match Python API
  } else {
    result.push_back(dycen); // duplicate to match Python API
    result.push_back(dycen); // duplicate to match Python API
  }

  return result;
}

// C++ implementation of the Gaussian profile function
py::array_t<double> gaussian(py::array_t<double> x, py::array_t<double> y,
                             double sky, double height, double xcen,
                             double ycen, double fwhm, int xbin, int ybin,
                             int ndiv) {

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  check_matching_shapes({&x_info, &y_info});

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output array with same shape as input
  py::array_t<double> result = py::array_t<double>(x_info.shape);
  py::buffer_info result_info = result.request();
  double *RESTRICT result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Gaussian profile parameter
  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  if (ndiv > 0) {
    // With sub-pixellation
    const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
    const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

// OpenMP parallelization of for loop with SIMD vectorization.
#ifdef _OPENMP
#pragma omp parallel for simd
#endif

    // Loop over all pixels
    for (std::size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double prof = 0.0;

      // Reuse precomputed sub-pixel offsets across all pixels.
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;
          prof += std::exp(-alpha * rsq);
        }
      }

      result_ptr[pixel_idx] = sky + height * inv_nadd * prof;
    }
  } else {
    // Fast calculation at pixel centers
    for (std::size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;
      result_ptr[i] = sky + height * std::exp(-alpha * rsq);
    }
  }

  return result;
}

// C++ implementation of the Gaussian profile derivatives
std::vector<py::array_t<double>>
dgaussian(py::array_t<double> x, py::array_t<double> y, double sky,
          double height, double xcen, double ycen, double fwhm, int xbin,
          int ybin, int ndiv, bool comp_dfwhm) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  check_matching_shapes({&x_info, &y_info});

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output arrays with same shape as input
  std::vector<py::array_t<double>> result;
  result.reserve(5);

  // Always need dsky, dheight, dxcen, dycen
  py::array_t<double> dsky = py::array_t<double>(x_info.shape);
  py::array_t<double> dheight = py::array_t<double>(x_info.shape);
  py::array_t<double> dxcen = py::array_t<double>(x_info.shape);
  py::array_t<double> dycen = py::array_t<double>(x_info.shape);

  // Optionally need dfwhm
  py::array_t<double> dfwhm;

  if (comp_dfwhm) {
    dfwhm = py::array_t<double>(x_info.shape);
  }

  // Get data pointers
  py::buffer_info dsky_info = dsky.request();
  py::buffer_info dheight_info = dheight.request();
  py::buffer_info dxcen_info = dxcen.request();
  py::buffer_info dycen_info = dycen.request();

  double *RESTRICT dsky_ptr = static_cast<double *>(dsky_info.ptr);
  double *RESTRICT dheight_ptr = static_cast<double *>(dheight_info.ptr);
  double *RESTRICT dxcen_ptr = static_cast<double *>(dxcen_info.ptr);
  double *RESTRICT dycen_ptr = static_cast<double *>(dycen_info.ptr);

  py::buffer_info dfwhm_info;
  double *RESTRICT dfwhm_ptr = nullptr;

  if (comp_dfwhm) {
    dfwhm_info = dfwhm.request();
    dfwhm_ptr = static_cast<double *>(dfwhm_info.ptr);
  }

  // Calculate Gaussian profile parameter
  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  double two_alpha_height = 2.0 * alpha * height;
  double dfwhm_coeff = two_alpha_height / fwhm;

  std::size_t n_pixels = static_cast<std::size_t>(x_info.size);

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  // Initialize dsky to ones (derivative of sky is always 1)
  std::fill_n(dsky_ptr, n_pixels, 1.0);

  if (ndiv > 0) {
    // With sub-pixellation
    const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
    const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);
    double inv_nadd =
        1.0 / static_cast<double>(x_offsets.size() * y_offsets.size());

// OpenMP parallelization of for loop with SIMD vectorization.
#ifdef _OPENMP
#pragma omp parallel for simd
#endif

    // Loop over all pixels
    for (std::size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      // Use _sum suffix to avoid shadowing the outer py::array_t declarations.
      double dheight_sum = 0.0;
      double dxcen_sum = 0.0;
      double dycen_sum = 0.0;
      double dfwhm_sum = 0.0;

      // Reuse precomputed sub-pixel offsets across all pixels.
      for (double yoff : y_offsets) {
        for (double xoff : x_offsets) {
          double dx = x_val + xoff - xcen;
          double dy = y_val + yoff - ycen;
          double rsq = dx * dx + dy * dy;

          // Gaussian value
          double dh = std::exp(-alpha * rsq);
          dheight_sum += dh;
          dxcen_sum += two_alpha_height * dh * dx;
          dycen_sum += two_alpha_height * dh * dy;

          if (comp_dfwhm) {
            dfwhm_sum += dfwhm_coeff * dh * rsq;
          }
        }
      }

      dheight_ptr[pixel_idx] = dheight_sum * inv_nadd;
      dxcen_ptr[pixel_idx] = dxcen_sum * inv_nadd;
      dycen_ptr[pixel_idx] = dycen_sum * inv_nadd;

      if (comp_dfwhm) {
        dfwhm_ptr[pixel_idx] = dfwhm_sum * inv_nadd;
      }
    }
  } else {
    // Fast calculation at pixel centers
    for (std::size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;

      // Gaussian value
      double dh = std::exp(-alpha * rsq);
      dheight_ptr[i] = dh;
      dxcen_ptr[i] = two_alpha_height * dh * dx;
      dycen_ptr[i] = two_alpha_height * dh * dy;

      if (comp_dfwhm) {
        dfwhm_ptr[i] = dfwhm_coeff * dh * rsq;
      }
    }
  }

  // Add arrays to result vector
  result.push_back(dsky);
  result.push_back(dheight);
  result.push_back(dxcen);
  result.push_back(dycen);

  if (comp_dfwhm) {
    result.push_back(dfwhm);
  } else {
    result.push_back(dycen); // duplicate to match Python API
  }

  return result;
}

namespace {

// Solves the weighted linear least squares problem for the fast Gaussian
// fitting (FGF) method of Wan et al. 2018, Sensors 18, 2836. For a symmetric
// Gaussian, ln(I) = ln(A) - m*((x-xc)**2 + (y-yc)**2), so that
//
//   I*(ln(I) + m*(x**2 + y**2) + p*x + q*y + k) = 0
//
// with p = -2*m*xc, q = -2*m*yc and k = m*(xc**2 + yc**2) - ln(A), which is
// linear in the unknowns. The unknowns are (m, p, q, k) or just (p, q, k) if
// m is held fixed (npar = 4 or 3). Each equation is divided by sigma. The
// normal matrix is returned in `cov` after inversion. Returns false if the
// system is singular.
//
// Arguments: xs, ys = pixel positions relative to the approximate centre;
// is = sky-subtracted counts (must be > 0 for the selected pixels); ss =
// per-pixel sigma; sel = indices of the pixels to use; fix_m / m_fixed = hold
// m = 4 ln2 / FWHM**2 at m_fixed instead of solving for it. On return theta
// holds the solution in the order above and cov the (npar x npar) inverse
// normal matrix, stored with a row stride of 4.
bool fgf_solve(const std::vector<double> &xs, const std::vector<double> &ys,
               const std::vector<double> &is, const std::vector<double> &ss,
               const std::vector<std::size_t> &sel, bool fix_m, double m_fixed,
               double *theta, double *cov) {
  const int npar = fix_m ? 3 : 4;
  // column offset of p in each row (column 0 is m when it is free)
  const int off = fix_m ? 0 : 1;
  if (sel.size() < static_cast<std::size_t>(npar)) {
    return false;
  }

  // Accumulate the normal equations A^T A theta = A^T b
  double ata[16] = {0.0};
  double atb[4] = {0.0};
  for (std::size_t idx : sel) {
    const double x = xs[idx], y = ys[idx], val = is[idx];
    // The I factor in the linearised equation (and 1/sigma) is applied as a
    // weight, which keeps the noise of ln(I) from dominating the fit
    const double w = val / ss[idx];
    double row[4];
    if (!fix_m) {
      row[0] = w * (x * x + y * y);
    }
    row[off] = w * x;
    row[off + 1] = w * y;
    row[off + 2] = w;
    double rhs = -w * std::log(val);
    if (fix_m) {
      // move the known m*(x**2 + y**2) term to the right-hand side
      rhs -= w * m_fixed * (x * x + y * y);
    }
    for (int i = 0; i < npar; ++i) {
      atb[i] += row[i] * rhs;
      for (int j = 0; j < npar; ++j) {
        ata[i * 4 + j] += row[i] * row[j];
      }
    }
  }

  // Gauss-Jordan elimination on [ata | identity] to get solution and inverse
  double aug[4][8] = {{0.0}};
  for (int i = 0; i < npar; ++i) {
    for (int j = 0; j < npar; ++j) {
      aug[i][j] = ata[i * 4 + j];
    }
    aug[i][npar + i] = 1.0;
  }
  for (int c = 0; c < npar; ++c) {
    // partial pivoting: pick the largest element in this column
    int piv = c;
    for (int r = c + 1; r < npar; ++r) {
      if (std::fabs(aug[r][c]) > std::fabs(aug[piv][c])) {
        piv = r;
      }
    }
    if (!(std::fabs(aug[piv][c]) > 1e-300)) {
      return false;
    }
    if (piv != c) {
      for (int j = 0; j < 2 * npar; ++j) {
        std::swap(aug[piv][j], aug[c][j]);
      }
    }
    const double d = aug[c][c];
    // normalise the pivot row, then eliminate this column from other rows
    for (int j = 0; j < 2 * npar; ++j) {
      aug[c][j] /= d;
    }
    for (int r = 0; r < npar; ++r) {
      if (r != c) {
        const double f = aug[r][c];
        for (int j = 0; j < 2 * npar; ++j) {
          aug[r][j] -= f * aug[c][j];
        }
      }
    }
  }
  // The right half of aug is now the inverse of A^T A (the covariance);
  // theta = inverse * A^T b
  for (int i = 0; i < npar; ++i) {
    double t = 0.0;
    for (int j = 0; j < npar; ++j) {
      t += aug[i][npar + j] * atb[j];
      cov[i * 4 + j] = aug[i][npar + j];
    }
    theta[i] = t;
  }
  return true;
}

} // namespace

namespace {

// Result of one two-step closed-form fit at a fixed sky level.
struct FgfResult {
  bool ok = false;
  double height = 0.0, xc = 0.0, yc = 0.0, fwhm = 0.0;
  double cov[16] = {0.0}; // covariance of (height, xc, yc, fwhm)
  double chisq = 0.0;     // chi**2 of the fit in the data domain, all pixels
};

// Pixel data shared by all fits, relative to the approximate centre.
struct FgfData {
  std::vector<double> xs, ys, raw, sg; // all OK pixels
  bool fix_fwhm = false;
  double m_fixed = 0.0;
  double snr_min = 5.0;
  bool direct_snr = true; // select pixels using data/sigma, single pass
};

// Two-step FGF of Wan et al. 2018 at sky level `sky`. Step 1 fits the five
// brightest pixels to estimate each pixel's SNR; step 2 fits all pixels with
// SNR > snr_min. If d.direct_snr, the SNR is instead taken directly as
// (data - sky)/sigma, so pixels are selected with no preliminary fit and a
// single closed-form solve is made.
FgfResult fgf_at_sky(const FgfData &d, double sky) {
  FgfResult res;
  const std::size_t n = d.xs.size();
  const int npar = d.fix_fwhm ? 3 : 4;

  // Sky-subtracted counts; only positive pixels have a defined logarithm
  std::vector<double> is(n);
  std::vector<std::size_t> pos;
  for (std::size_t j = 0; j < n; ++j) {
    is[j] = d.raw[j] - sky;
    if (is[j] > 0.0) {
      pos.push_back(j);
    }
  }
  // need more pixels than parameters for the fit to mean anything
  if (pos.size() < static_cast<std::size_t>(npar) + 1) {
    return res;
  }

  double theta[4], cov[16];
  // Fallback selection: the (up to) 5 brightest pixels
  std::vector<std::size_t> sel = pos;
  const std::size_t nbright = std::min<std::size_t>(5, sel.size());
  std::partial_sort(sel.begin(), sel.begin() + nbright, sel.end(),
                    [&](std::size_t a, std::size_t b) { return is[a] > is[b]; });
  sel.resize(nbright);

  // Converts the solved linear parameters back to (m, p, q, k)
  auto unpack = [&](double &m, double &p, double &q, double &k) {
    if (d.fix_fwhm) {
      m = d.m_fixed;
      p = theta[0];
      q = theta[1];
      k = theta[2];
    } else {
      m = theta[0];
      p = theta[1];
      q = theta[2];
      k = theta[3];
    }
  };

  double m, p, q, k;
  std::vector<std::size_t> sel2;
  if (d.direct_snr) {
    // keep pixels whose own SNR exceeds the threshold
    for (std::size_t j : pos) {
      if (is[j] > d.snr_min * d.sg[j]) {
        sel2.push_back(j);
      }
    }
  } else {
    // Step 1: rough fit to the brightest pixels
    if (!fgf_solve(d.xs, d.ys, is, d.sg, sel, d.fix_fwhm, d.m_fixed, theta,
                   cov)) {
      return res;
    }
    unpack(m, p, q, k);
    // m <= 0 would be a rising, not a falling, profile
    if (!(m > 0.0)) {
      return res;
    }

    // Step 2 selection: model signal s from the rough fit; keep pixels where
    // the signal exceeds snr_min times the residual (data - model)
    const double lnA1 = (p * p + q * q) / (4.0 * m) - k;
    const double xc1 = -p / (2.0 * m), yc1 = -q / (2.0 * m);
    for (std::size_t j : pos) {
      const double dx = d.xs[j] - xc1, dy = d.ys[j] - yc1;
      const double s = std::exp(lnA1 - m * (dx * dx + dy * dy));
      if (std::fabs(s) > d.snr_min * std::fabs(is[j] - s)) {
        sel2.push_back(j);
      }
    }
  }
  // too few pixels survived the SNR cut: fall back to the brightest pixels
  if (sel2.size() < static_cast<std::size_t>(npar) + 1) {
    sel2 = sel;
  }
  // Final closed-form solve using the selected pixels
  if (!fgf_solve(d.xs, d.ys, is, d.sg, sel2, d.fix_fwhm, d.m_fixed, theta,
                 cov)) {
    return res;
  }
  unpack(m, p, q, k);
  if (!(m > 0.0)) {
    return res;
  }

  // Convert (m, p, q, k) to physical parameters: centre, amplitude (peak
  // height above sky) and FWHM
  const double xc = -p / (2.0 * m), yc = -q / (2.0 * m);
  const double A = std::exp((p * p + q * q) / (4.0 * m) - k);
  const double fw = std::sqrt(4.0 * std::log(2.0) / m);
  if (!std::isfinite(A) || !std::isfinite(fw)) {
    return res;
  }

  // Jacobian of (A, xc, yc, fwhm) w.r.t. the solved parameters
  const int off = d.fix_fwhm ? 0 : 1;
  double jac[4][4] = {{0.0}};
  if (!d.fix_fwhm) {
    jac[0][0] = -A * (p * p + q * q) / (4.0 * m * m);
    jac[1][0] = p / (2.0 * m * m);
    jac[2][0] = q / (2.0 * m * m);
    jac[3][0] = -fw / (2.0 * m);
  }
  jac[0][off] = A * p / (2.0 * m);
  jac[0][off + 1] = A * q / (2.0 * m);
  jac[0][off + 2] = -A;
  jac[1][off] = -1.0 / (2.0 * m);
  jac[2][off + 1] = -1.0 / (2.0 * m);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      double s = 0.0;
      for (int i = 0; i < npar; ++i) {
        for (int j = 0; j < npar; ++j) {
          s += jac[a][i] * cov[i * 4 + j] * jac[b][j];
        }
      }
      res.cov[a * 4 + b] = s;
    }
  }

  // chi**2 in the data domain over all pixels, used to fit the sky
  double chisq = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    const double dx = d.xs[j] - xc, dy = d.ys[j] - yc;
    const double r = (is[j] - A * std::exp(-m * (dx * dx + dy * dy))) / d.sg[j];
    chisq += r * r;
  }
  res.ok = true;
  res.height = A;
  res.xc = xc;
  res.yc = yc;
  res.fwhm = fw;
  res.chisq = chisq;
  return res;
}

} // namespace

// Fast Gaussian fitting of a symmetric 2D Gaussian plus constant after Wan et
// al. 2018, Sensors 18, 2836, extended to fit the sky. For a given sky the fit
// is closed-form (see fgf_at_sky); if the sky is free it is found with a
// golden-section search on the resulting chi**2 over
// [sky - sky_range, sky + sky_range], stopping when the bracket is < sky_tol.
// The sky variance comes from the curvature of the profile chi**2.
//
// Returns 4 + 16 + 2 values: fitted (height, xcen, ycen, fwhm), their 4x4
// row-major covariance (unscaled apart from the input sigmas; ignores the
// correlation with the sky), then the sky and its variance (-1 if the sky was
// fixed). If fix_fwhm, the fwhm row/column of the covariance is zero.
//
// Arguments
// ---------
// x, y        : pixel coordinates of every pixel in the fit window (any
//               shape, but x, y, data and sigma must all match). Flattened
//               in C order, so indices in ok_indices refer to the flat array.
// data        : debiassed counts in each pixel, including the sky.
// sigma       : 1-sigma uncertainty of each pixel, in counts. Used to weight
//               the fit and to select pixels by SNR.
// ok_indices  : 1D int64 array of flat indices of the pixels to use. Allows
//               masked/bad/saturated pixels to be left out.
// sky         : if fix_sky, the sky level to subtract; otherwise the centre
//               of the range searched for the best sky.
// fix_sky     : true to hold the sky at `sky`; false to fit it.
// sky_range   : half-width of the sky search interval (counts); only used if
//               fix_sky is false. Must be > 0 in that case.
// sky_tol     : stop the sky search when the bracket is narrower than this
//               (counts); only used if fix_sky is false. Must be > 0 in that
//               case.
// xcen, ycen  : approximate centre of the profile, in the same coordinates as
//               x and y. Coordinates are measured relative to this point
//               internally for numerical stability. The fit does not need an
//               accurate position, but it should be reasonably close.
// fwhm        : FWHM in pixels. Held fixed if fix_fwhm, otherwise only
//               checked to be positive (the fitted value is independent of
//               it).
// fix_fwhm    : true to hold the FWHM at `fwhm` and fit only the height and
//               centre; false to fit the FWHM too.
// snr_min     : minimum signal-to-noise for a pixel to be included in the
//               final fit. Faint pixels have the most noisy logarithms so
//               are excluded.
// direct_snr  : if true, SNR = (data - sky)/sigma, which needs just a single
//               closed-form solve. If false, the SNR is estimated from a
//               preliminary fit to the 5 brightest pixels as in Wan et al.
//               (two solves).
//
// Throws std::runtime_error if the inputs are invalid or the fit fails
// (e.g. too few usable pixels, or a non-positive Gaussian curvature).
std::vector<double> fast_gaussian_fit(
    py::array_t<double> x, py::array_t<double> y, py::array_t<double> data,
    py::array_t<double> sigma, py::array_t<std::int64_t> ok_indices,
    double sky, bool fix_sky, double sky_range, double sky_tol, double xcen,
    double ycen, double fwhm, bool fix_fwhm, double snr_min,
    bool direct_snr) {
  py::buffer_info xi = x.request(), yi = y.request(), di = data.request(),
                  si = sigma.request(), oi = ok_indices.request();
  check_matching_shapes({&xi, &yi, &di, &si});
  check_ok_indices_1d(oi);
  const std::size_t npix = static_cast<std::size_t>(xi.size);
  const std::size_t nok = static_cast<std::size_t>(oi.shape[0]);
  const auto *ok_ptr = static_cast<const std::int64_t *>(oi.ptr);
  check_ok_indices_range(ok_ptr, nok, npix);
  if (!(fwhm > 0.0)) {
    throw std::runtime_error("fwhm must be positive");
  }
  if (!fix_sky && !(sky_range > 0.0 && sky_tol > 0.0)) {
    throw std::runtime_error("sky_range and sky_tol must be positive");
  }

  const double *xp = static_cast<const double *>(xi.ptr);
  const double *yp = static_cast<const double *>(yi.ptr);
  const double *dp = static_cast<const double *>(di.ptr);
  const double *sp = static_cast<const double *>(si.ptr);

  // Gather the usable pixels into contiguous arrays, with positions relative
  // to the approximate centre, and convert the FWHM to the Gaussian exponent
  // m = 4 ln2 / FWHM**2 used when it is held fixed.
  FgfData d;
  d.fix_fwhm = fix_fwhm;
  d.m_fixed = 4.0 * std::log(2.0) / (fwhm * fwhm);
  d.snr_min = snr_min;
  d.direct_snr = direct_snr;
  for (std::size_t n = 0; n < nok; ++n) {
    const std::size_t idx = static_cast<std::size_t>(ok_ptr[n]);
    d.xs.push_back(xp[idx] - xcen);
    d.ys.push_back(yp[idx] - ycen);
    d.raw.push_back(dp[idx]);
    d.sg.push_back(std::fabs(sp[idx]));
  }

  double best = sky;
  double sky_var = -1.0;
  FgfResult res;
  if (fix_sky) {
    // Sky known: a single closed-form fit does everything
    res = fgf_at_sky(d, sky);
  } else {
    // Sky free: minimise the chi**2 of the closed-form fit over the sky
    // large finite sentinel: infinities are unsafe under -ffast-math
    const double kInfeasible = 1.0e300;
    auto f = [&](double c) {
      FgfResult r = fgf_at_sky(d, c);
      return r.ok ? r.chisq : kInfeasible;
    };
    // Coarse scan to bracket the minimum among feasible skies, because the
    // pixel selection changes with the sky and f can be infinite or jumpy.
    const int ngrid = 21;
    const double step = 2.0 * sky_range / (ngrid - 1);
    int ibest = -1;
    double fbest = kInfeasible;
    for (int i = 0; i < ngrid; ++i) {
      const double fi = f(sky - sky_range + i * step);
      if (fi < fbest) {
        fbest = fi;
        ibest = i;
      }
    }
    if (ibest < 0) {
      throw std::runtime_error("fast Gaussian fit failed");
    }

    // Golden-section refinement within the grid cells either side of the
    // best grid point. gr is the inverse golden ratio; c1 < c2 are the two
    // interior probe points of the bracket [a, b].
    const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
    double a = sky - sky_range + std::max(ibest - 1, 0) * step;
    double b = sky - sky_range + std::min(ibest + 1, ngrid - 1) * step;
    double c1 = b - gr * (b - a), c2 = a + gr * (b - a);
    double f1 = f(c1), f2 = f(c2);
    const double xbest = sky - sky_range + ibest * step;
    for (int it = 0; it < 100 && (b - a) > sky_tol; ++it) {
      // if both are infeasible, move towards the best grid point
      const bool left = (f1 >= kInfeasible && f2 >= kInfeasible) ? (xbest < c2) : f1 < f2;
      if (left) {
        // minimum lies in [a, c2]: shrink from the right, reusing c1
        b = c2;
        c2 = c1;
        f2 = f1;
        c1 = b - gr * (b - a);
        f1 = f(c1);
      } else {
        // minimum lies in [c1, b]: shrink from the left, reusing c2
        a = c1;
        c1 = c2;
        f1 = f2;
        c2 = a + gr * (b - a);
        f2 = f(c2);
      }
    }
    best = 0.5 * (a + b);
    res = fgf_at_sky(d, best);

    // Curvature of the profile chi**2 for the sky variance, from a central
    // second difference. For chi**2 = chi2_min + (sky - best)**2/var the
    // curvature is 2/var, hence var = 2/curv. Left at -1 if not measurable.
    const double h = std::max(sky_tol, 0.05 * sky_range);
    const double fp = f(best + h), fm = f(best - h);
    const double curv = (fp - 2.0 * res.chisq + fm) / (h * h);
    if (fp < kInfeasible && fm < kInfeasible && curv > 0.0) {
      sky_var = 2.0 / curv;
    }
  }
  if (!res.ok) {
    throw std::runtime_error("fast Gaussian fit failed");
  }

  // Pack the output: [0-3] height, xcen, ycen, fwhm; [4-19] covariance;
  // [20] sky; [21] sky variance. Positions are shifted back to input coords.
  std::vector<double> result(22, 0.0);
  result[0] = res.height;
  result[1] = xcen + res.xc;
  result[2] = ycen + res.yc;
  result[3] = res.fwhm;
  std::copy(res.cov, res.cov + 16, result.begin() + 4);
  result[20] = best;
  result[21] = fix_sky ? -1.0 : sky_var;
  return result;
}

PYBIND11_MODULE(_fitting_cpp, m) {
  m.doc() = "C++ implementation of profile fitting functions";

  m.def("fast_gaussian_fit", &fast_gaussian_fit,
        "Closed-form fast Gaussian fit (Wan et al. 2018, Sensors 18, 2836)",
        py::arg("x"), py::arg("y"), py::arg("data"), py::arg("sigma"),
        py::arg("ok_indices"), py::arg("sky"), py::arg("fix_sky"),
        py::arg("sky_range"), py::arg("sky_tol"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("fix_fwhm"),
        py::arg("snr_min") = 5.0, py::arg("direct_snr") = true);

  m.def("moffat", &moffat, "C++ implementation of Moffat profile", py::arg("x"),
        py::arg("y"), py::arg("sky"), py::arg("height"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("beta"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dmoffat", &dmoffat, "C++ implementation of Moffat profile derivatives",
        py::arg("x"), py::arg("y"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("beta"),
        py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"),
        py::arg("comp_dfwhm"), py::arg("comp_dbeta"));

  m.def("moffat_resid", &moffat_resid,
        "C++ implementation of Moffat residuals at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("data"), py::arg("sigma"),
        py::arg("ok_indices"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("beta"),
        py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"));

  m.def("dmoffat_jac", &dmoffat_jac,
        "C++ implementation of Moffat residual Jacobian at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("sigma"), py::arg("ok_indices"),
        py::arg("sky"), py::arg("height"), py::arg("xcen"), py::arg("ycen"),
        py::arg("fwhm"), py::arg("beta"), py::arg("xbin"), py::arg("ybin"),
        py::arg("ndiv"), py::arg("comp_dfwhm"), py::arg("comp_dbeta"),
        py::arg("inds"));

  m.def("gaussian", &gaussian, "C++ implementation of Gaussian profile",
        py::arg("x"), py::arg("y"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dgaussian", &dgaussian,
        "C++ implementation of Gaussian profile derivatives", py::arg("x"),
        py::arg("y"), py::arg("sky"), py::arg("height"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"), py::arg("ybin"),
        py::arg("ndiv"), py::arg("comp_dfwhm"));

  m.def("gaussian_resid", &gaussian_resid,
        "C++ implementation of Gaussian residuals at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("data"), py::arg("sigma"),
        py::arg("ok_indices"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dgaussian_jac", &dgaussian_jac,
        "C++ implementation of Gaussian residual Jacobian at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("sigma"), py::arg("ok_indices"),
        py::arg("sky"), py::arg("height"), py::arg("xcen"), py::arg("ycen"),
        py::arg("fwhm"), py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"),
        py::arg("comp_dfwhm"), py::arg("inds"));
}
