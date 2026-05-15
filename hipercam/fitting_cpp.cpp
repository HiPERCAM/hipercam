#include <algorithm>
#include <cmath>
#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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

// Helper function to generate sub-pixel offsets for binning.
// This creates a vector of offsets for each bin and sub-bin, centered around
// zero, to be used in the sub-pixellation loops in the profile evaluation and
// derivatives.
inline std::vector<double> make_subpixel_offsets(int bin, int ndiv) {
  std::vector<double> offsets;
  if (ndiv <= 0) {
    return offsets;
  }
  offsets.reserve(static_cast<size_t>(bin * ndiv));
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
inline double calc_moffat_alpha(double fwhm, double beta) {
  double tbeta = std::max(0.01, beta);
  return 4.0 * (std::pow(2.0, 1.0 / tbeta) - 1.0) / (fwhm * fwhm);
}

// Evaluate one Moffat model value at a single pixel coordinate.
// This is used by the selected-pixel residual path to avoid building
// full 2D model arrays when only masked pixels are needed.
inline double moffat_value_at(double x_val, double y_val, double height,
                              double xcen, double ycen, double alpha,
                              double tbeta,
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
inline double gaussian_value_at(double x_val, double y_val, double height,
                                double xcen, double ycen, double alpha,
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
// The outputs follow the same normalization as dmoffat_cpp so the fit path
// stays numerically equivalent.
inline void
moffat_derivs_at(double x_val, double y_val, double height, double xcen,
                 double ycen, double alpha, double tbeta, double dfwhm_coeff,
                 double dbeta_coeff, const std::vector<double> &x_offsets,
                 const std::vector<double> &y_offsets, bool comp_dfwhm,
                 bool comp_dbeta, double &dheight, double &dxcen, double &dycen,
                 double &dfwhm, double &dbeta) {
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
// The outputs follow the same normalization as dgaussian_cpp so the fit path
// stays numerically equivalent.
inline void gaussian_derivs_at(double x_val, double y_val, double height,
                               double xcen, double ycen, double alpha,
                               double two_alpha_height, double dfwhm_coeff,
                               const std::vector<double> &x_offsets,
                               const std::vector<double> &y_offsets,
                               bool comp_dfwhm, double &dheight, double &dxcen,
                               double &dycen, double &dfwhm) {
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

py::array_t<double>
moffat_resid_cpp(py::array_t<double> x, py::array_t<double> y,
                 py::array_t<double> data, py::array_t<double> sigma,
                 py::array_t<std::int64_t> ok_indices, double sky,
                 double height, double xcen, double ycen, double fwhm,
                 double beta, int xbin, int ybin, int ndiv) {

  // Compute residuals only at valid indices provided by Python. This
  // bypasses full-frame residual assembly and boolean masking in Python.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info data_info = data.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 || data_info.ndim != 2 ||
      sigma_info.ndim != 2 || ok_info.ndim != 1 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1] ||
      x_info.shape[0] != data_info.shape[0] ||
      x_info.shape[1] != data_info.shape[1] ||
      x_info.shape[0] != sigma_info.shape[0] ||
      x_info.shape[1] != sigma_info.shape[1]) {
    throw std::runtime_error(
        "Input arrays have invalid dimensions or mismatched shapes");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *data_ptr = static_cast<double *>(data_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];
  size_t n_ok = ok_info.shape[0];

  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(static_cast<py::ssize_t>(n_ok));
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];
    if (idx < 0 || static_cast<size_t>(idx) >= n_pixels) {
      throw std::runtime_error("ok_indices contains out-of-range values");
    }
  }

  // Release the GIL now the Python operations are complete.
  // pybind11 buffer operations (request, array creation) require the GIL,
  // so we can only release it after extracting all pointers and dimensions.
  py::gil_scoped_release release;

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double model =
        sky + moffat_value_at(x_ptr[idx], y_ptr[idx], height, xcen, ycen, alpha,
                              tbeta, x_offsets, y_offsets);
    result_ptr[i] = (data_ptr[idx] - model) / sigma_ptr[idx];
  }

  return result;
}

py::array_t<double> dmoffat_jac_cpp(
    py::array_t<double> x, py::array_t<double> y, py::array_t<double> sigma,
    py::array_t<std::int64_t> ok_indices, double sky, double height,
    double xcen, double ycen, double fwhm, double beta, int xbin, int ybin,
    int ndiv, bool comp_dfwhm, bool comp_dbeta, const std::vector<int> &inds) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Build Jacobian rows directly for selected pixels and selected parameter
  // columns (`inds`), matching the same derivative ordering used by Python.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 || sigma_info.ndim != 2 ||
      ok_info.ndim != 1 || x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1] ||
      x_info.shape[0] != sigma_info.shape[0] ||
      x_info.shape[1] != sigma_info.shape[1]) {
    throw std::runtime_error(
        "Input arrays have invalid dimensions or mismatched shapes");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];
  size_t n_ok = ok_info.shape[0];
  size_t n_par = inds.size();

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

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];
    if (idx < 0 || static_cast<size_t>(idx) >= n_pixels) {
      throw std::runtime_error("ok_indices contains out-of-range values");
    }
  }

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (size_t i = 0; i < n_ok; ++i) {
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

    for (size_t j = 0; j < n_par; ++j) {
      int ind = inds[j];
      if (ind < 0 || ind > 5) {
        throw std::runtime_error("inds contains out-of-range derivative index");
      }
      result_ptr[i * n_par + j] = derivs[ind] * inv_sigma;
    }
  }

  return result;
}

py::array_t<double>
gaussian_resid_cpp(py::array_t<double> x, py::array_t<double> y,
                   py::array_t<double> data, py::array_t<double> sigma,
                   py::array_t<std::int64_t> ok_indices, double sky,
                   double height, double xcen, double ycen, double fwhm,
                   int xbin, int ybin, int ndiv) {

  // Gaussian equivalent of moffat_resid_cpp: selected-pixel residuals only.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info data_info = data.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 || data_info.ndim != 2 ||
      sigma_info.ndim != 2 || ok_info.ndim != 1 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1] ||
      x_info.shape[0] != data_info.shape[0] ||
      x_info.shape[1] != data_info.shape[1] ||
      x_info.shape[0] != sigma_info.shape[0] ||
      x_info.shape[1] != sigma_info.shape[1]) {
    throw std::runtime_error(
        "Input arrays have invalid dimensions or mismatched shapes");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *data_ptr = static_cast<double *>(data_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];
  size_t n_ok = ok_info.shape[0];

  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(static_cast<py::ssize_t>(n_ok));
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];
    if (idx < 0 || static_cast<size_t>(idx) >= n_pixels) {
      throw std::runtime_error("ok_indices contains out-of-range values");
    }
  }

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];

    double model = sky + gaussian_value_at(x_ptr[idx], y_ptr[idx], height, xcen,
                                           ycen, alpha, x_offsets, y_offsets);
    result_ptr[i] = (data_ptr[idx] - model) / sigma_ptr[idx];
  }

  return result;
}

py::array_t<double> dgaussian_jac_cpp(
    py::array_t<double> x, py::array_t<double> y, py::array_t<double> sigma,
    py::array_t<std::int64_t> ok_indices, double sky, double height,
    double xcen, double ycen, double fwhm, int xbin, int ybin, int ndiv,
    bool comp_dfwhm, const std::vector<int> &inds) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Selected-pixel Jacobian assembly for Gaussian fits.

  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();
  py::buffer_info sigma_info = sigma.request();
  py::buffer_info ok_info = ok_indices.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 || sigma_info.ndim != 2 ||
      ok_info.ndim != 1 || x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1] ||
      x_info.shape[0] != sigma_info.shape[0] ||
      x_info.shape[1] != sigma_info.shape[1]) {
    throw std::runtime_error(
        "Input arrays have invalid dimensions or mismatched shapes");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);
  const double *sigma_ptr = static_cast<double *>(sigma_info.ptr);
  const std::int64_t *ok_ptr = static_cast<std::int64_t *>(ok_info.ptr);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];
  size_t n_ok = ok_info.shape[0];
  size_t n_par = inds.size();

  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  double two_alpha_height = 2.0 * alpha * height;
  double dfwhm_coeff = two_alpha_height / fwhm;
  const std::vector<double> x_offsets = make_subpixel_offsets(xbin, ndiv);
  const std::vector<double> y_offsets = make_subpixel_offsets(ybin, ndiv);

  py::array_t<double> result(
      {static_cast<py::ssize_t>(n_ok), static_cast<py::ssize_t>(n_par)});
  py::buffer_info result_info = result.request();
  double *result_ptr = static_cast<double *>(result_info.ptr);

  for (size_t i = 0; i < n_ok; ++i) {
    std::int64_t idx = ok_ptr[i];
    if (idx < 0 || static_cast<size_t>(idx) >= n_pixels) {
      throw std::runtime_error("ok_indices contains out-of-range values");
    }
  }

  // Release the GIL now the Python operations are complete.
  py::gil_scoped_release release;

  for (size_t i = 0; i < n_ok; ++i) {
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

    for (size_t j = 0; j < n_par; ++j) {
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
py::array_t<double> moffat_cpp(py::array_t<double> x, py::array_t<double> y,
                               double sky, double height, double xcen,
                               double ycen, double fwhm, double beta, int xbin,
                               int ybin, int ndiv) {

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1]) {
    throw std::runtime_error("Input arrays must be 2D and have the same shape");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output array with same shape as input
  py::array_t<double> result = py::array_t<double>(x_info.shape);
  py::buffer_info result_info = result.request();
  double *RESTRICT result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Moffat profile parameters
  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

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
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
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
    for (size_t i = 0; i < n_pixels; ++i) {
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
dmoffat_cpp(py::array_t<double> x, py::array_t<double> y, double sky,
            double height, double xcen, double ycen, double fwhm, double beta,
            int xbin, int ybin, int ndiv, bool comp_dfwhm, bool comp_dbeta) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1]) {
    throw std::runtime_error("Input arrays must be 2D and have the same shape");
  }

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

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

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
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
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
    for (size_t i = 0; i < n_pixels; ++i) {
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
py::array_t<double> gaussian_cpp(py::array_t<double> x, py::array_t<double> y,
                                 double sky, double height, double xcen,
                                 double ycen, double fwhm, int xbin, int ybin,
                                 int ndiv) {

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1]) {
    throw std::runtime_error("Input arrays must be 2D and have the same shape");
  }

  const double *x_ptr = static_cast<double *>(x_info.ptr);
  const double *y_ptr = static_cast<double *>(y_info.ptr);

  // Create output array with same shape as input
  py::array_t<double> result = py::array_t<double>(x_info.shape);
  py::buffer_info result_info = result.request();
  double *RESTRICT result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Gaussian profile parameter
  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

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
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
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
    for (size_t i = 0; i < n_pixels; ++i) {
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
dgaussian_cpp(py::array_t<double> x, py::array_t<double> y, double sky,
              double height, double xcen, double ycen, double fwhm, int xbin,
              int ybin, int ndiv, bool comp_dfwhm) {

  // sky affects residuals but not derivatives; suppress unused-param warning
  (void)sky;

  // Get input array dimensions and data pointers
  py::buffer_info x_info = x.request();
  py::buffer_info y_info = y.request();

  if (x_info.ndim != 2 || y_info.ndim != 2 ||
      x_info.shape[0] != y_info.shape[0] ||
      x_info.shape[1] != y_info.shape[1]) {
    throw std::runtime_error("Input arrays must be 2D and have the same shape");
  }

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

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

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
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
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
    for (size_t i = 0; i < n_pixels; ++i) {
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

PYBIND11_MODULE(fitting_cpp, m) {
  m.doc() = "C++ implementation of profile fitting functions";

  m.def("moffat", &moffat_cpp, "C++ implementation of Moffat profile",
        py::arg("x"), py::arg("y"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("beta"),
        py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"));

  m.def("dmoffat", &dmoffat_cpp,
        "C++ implementation of Moffat profile derivatives", py::arg("x"),
        py::arg("y"), py::arg("sky"), py::arg("height"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("beta"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"), py::arg("comp_dfwhm"),
        py::arg("comp_dbeta"));

  m.def("moffat_resid", &moffat_resid_cpp,
        "C++ implementation of Moffat residuals at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("data"), py::arg("sigma"),
        py::arg("ok_indices"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("beta"),
        py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"));

  m.def("dmoffat_jac", &dmoffat_jac_cpp,
        "C++ implementation of Moffat residual Jacobian at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("sigma"), py::arg("ok_indices"),
        py::arg("sky"), py::arg("height"), py::arg("xcen"), py::arg("ycen"),
        py::arg("fwhm"), py::arg("beta"), py::arg("xbin"), py::arg("ybin"),
        py::arg("ndiv"), py::arg("comp_dfwhm"), py::arg("comp_dbeta"),
        py::arg("inds"));

  m.def("gaussian", &gaussian_cpp, "C++ implementation of Gaussian profile",
        py::arg("x"), py::arg("y"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dgaussian", &dgaussian_cpp,
        "C++ implementation of Gaussian profile derivatives", py::arg("x"),
        py::arg("y"), py::arg("sky"), py::arg("height"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"), py::arg("ybin"),
        py::arg("ndiv"), py::arg("comp_dfwhm"));

  m.def("gaussian_resid", &gaussian_resid_cpp,
        "C++ implementation of Gaussian residuals at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("data"), py::arg("sigma"),
        py::arg("ok_indices"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dgaussian_jac", &dgaussian_jac_cpp,
        "C++ implementation of Gaussian residual Jacobian at selected pixels",
        py::arg("x"), py::arg("y"), py::arg("sigma"), py::arg("ok_indices"),
        py::arg("sky"), py::arg("height"), py::arg("xcen"), py::arg("ycen"),
        py::arg("fwhm"), py::arg("xbin"), py::arg("ybin"), py::arg("ndiv"),
        py::arg("comp_dfwhm"), py::arg("inds"));
}
