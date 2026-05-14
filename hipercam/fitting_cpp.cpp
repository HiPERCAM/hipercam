#include <algorithm>
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

// Helper function to calculate the Moffat alpha parameter
inline double calc_moffat_alpha(double fwhm, double beta) {
  double tbeta = std::max(0.01, beta);
  return 4.0 * (std::pow(2.0, 1.0 / tbeta) - 1.0) / (fwhm * fwhm);
}

/**
 * C++ implementation of the Moffat profile function
 */
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
  double *result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Moffat profile parameters
  double tbeta = std::max(0.01, beta);
  double alpha = calc_moffat_alpha(fwhm, beta);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

  if (ndiv > 0) {
    // With sub-pixellation
    std::fill_n(result_ptr, n_pixels, 0.0);
    double norm = height / xbin / ybin / (ndiv * ndiv);
    double inv_ndiv = 1.0 / static_cast<double>(ndiv);

    // Mean offset within sub-pixels
    double soff = (ndiv - 1.0) / (2.0 * ndiv);

    // Loop over all pixels
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double prof = 0.0;

      // Loop over sub-pixels
      for (int iy = 0; iy < ybin; ++iy) {
        double yoff = iy - (ybin - 1) / 2.0 - soff;
        for (int ix = 0; ix < xbin; ++ix) {
          double xoff = ix - (xbin - 1) / 2.0 - soff;
          for (int isy = 0; isy < ndiv; ++isy) {
            double ysoff = yoff + isy * inv_ndiv;
            for (int isx = 0; isx < ndiv; ++isx) {
              double xsoff = xoff + isx * inv_ndiv;
              double dx = x_val + xsoff - xcen;
              double dy = y_val + ysoff - ycen;
              double rsq = dx * dx + dy * dy;
              prof += norm * std::pow(1.0 + alpha * rsq, -tbeta);
            }
          }
        }
      }

      result_ptr[pixel_idx] = sky + prof;
    }
  } else {
    // Fast calculation at pixel centers
    for (size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;
      result_ptr[i] = height * std::pow(1.0 + alpha * rsq, -tbeta) + sky;
    }
  }

  return result;
}

/**
 * C++ implementation of the Moffat profile derivatives
 */
std::vector<py::array_t<double>>
dmoffat_cpp(py::array_t<double> x, py::array_t<double> y, double sky,
            double height, double xcen, double ycen, double fwhm, double beta,
            int xbin, int ybin, int ndiv, bool comp_dfwhm, bool comp_dbeta) {

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

  double *dsky_ptr = static_cast<double *>(dsky_info.ptr);
  double *dheight_ptr = static_cast<double *>(dheight_info.ptr);
  double *dxcen_ptr = static_cast<double *>(dxcen_info.ptr);
  double *dycen_ptr = static_cast<double *>(dycen_info.ptr);

  py::buffer_info dfwhm_info;
  py::buffer_info dbeta_info;
  double *dfwhm_ptr = nullptr;
  double *dbeta_ptr = nullptr;

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

  // Initialize dsky to ones (derivative of sky is always 1)
  std::fill_n(dsky_ptr, n_pixels, 1.0);

  if (ndiv > 0) {
    // With sub-pixellation
    double inv_nadd = 1.0 / static_cast<double>(xbin * ybin * ndiv * ndiv);
    double inv_ndiv = 1.0 / static_cast<double>(ndiv);

    // Mean offset within sub-pixels
    double soff = (ndiv - 1.0) / (2.0 * ndiv);

    // Loop over all pixels
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double dheight = 0.0;
      double dxcen = 0.0;
      double dycen = 0.0;
      double dfwhm = 0.0;
      double dbeta = 0.0;

      // Loop over sub-pixels
      for (int iy = 0; iy < ybin; ++iy) {
        double yoff = iy - (ybin - 1) / 2.0 - soff;
        for (int ix = 0; ix < xbin; ++ix) {
          double xoff = ix - (xbin - 1) / 2.0 - soff;
          for (int isy = 0; isy < ndiv; ++isy) {
            double ysoff = yoff + isy * inv_ndiv;
            for (int isx = 0; isx < ndiv; ++isx) {
              double xsoff = xoff + isx * inv_ndiv;
              double dx = x_val + xsoff - xcen;
              double dy = y_val + ysoff - ycen;
              double rsq = dx * dx + dy * dy;

              double denom = 1.0 + alpha * rsq;
              double save1 = height * std::pow(denom, -tbeta - 1.0);
              double save2 = save1 * rsq;

              // Derivatives
              double dh = std::pow(denom, -tbeta);
              dheight += dh;
              dxcen += two_alpha_tbeta * dx * save1;
              dycen += two_alpha_tbeta * dy * save1;

              if (comp_dfwhm) {
                dfwhm += dfwhm_coeff * save2;
              }

              if (comp_dbeta) {
                double log_denom = std::log(denom);
                dbeta += (-log_denom * height * dh + dbeta_coeff * save2);
              }
            }
          }
        }
      }

      dheight_ptr[pixel_idx] = dheight * inv_nadd;
      dxcen_ptr[pixel_idx] = dxcen * inv_nadd;
      dycen_ptr[pixel_idx] = dycen * inv_nadd;

      if (comp_dfwhm) {
        dfwhm_ptr[pixel_idx] = dfwhm * inv_nadd;
      }

      if (comp_dbeta) {
        dbeta_ptr[pixel_idx] = dbeta * inv_nadd;
      }
    }
  } else {
    // Fast calculation at pixel centers
    for (size_t i = 0; i < n_pixels; ++i) {
      double dx = x_ptr[i] - xcen;
      double dy = y_ptr[i] - ycen;
      double rsq = dx * dx + dy * dy;

      double denom = 1.0 + alpha * rsq;
      double save1 = height * std::pow(denom, -tbeta - 1.0);
      double save2 = save1 * rsq;

      // Derivatives
      dheight_ptr[i] = std::pow(denom, -tbeta);
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

/**
 * C++ implementation of the Gaussian profile function
 */
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
  double *result_ptr = static_cast<double *>(result_info.ptr);

  // Calculate Gaussian profile parameter
  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

  if (ndiv > 0) {
    // With sub-pixellation
    std::fill_n(result_ptr, n_pixels, 0.0);
    double norm = height / xbin / ybin / (ndiv * ndiv);
    double inv_ndiv = 1.0 / static_cast<double>(ndiv);

    // Mean offset within sub-pixels
    double soff = (ndiv - 1.0) / (2.0 * ndiv);

    // Loop over all pixels
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double prof = 0.0;

      // Loop over sub-pixels
      for (int iy = 0; iy < ybin; ++iy) {
        double yoff = iy - (ybin - 1) / 2.0 - soff;
        for (int ix = 0; ix < xbin; ++ix) {
          double xoff = ix - (xbin - 1) / 2.0 - soff;
          for (int isy = 0; isy < ndiv; ++isy) {
            double ysoff = yoff + isy * inv_ndiv;
            for (int isx = 0; isx < ndiv; ++isx) {
              double xsoff = xoff + isx * inv_ndiv;
              double dx = x_val + xsoff - xcen;
              double dy = y_val + ysoff - ycen;
              double rsq = dx * dx + dy * dy;
              prof += std::exp(-alpha * rsq);
            }
          }
        }
      }

      result_ptr[pixel_idx] = sky + norm * prof;
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

/**
 * C++ implementation of the Gaussian profile derivatives
 */
std::vector<py::array_t<double>>
dgaussian_cpp(py::array_t<double> x, py::array_t<double> y, double sky,
              double height, double xcen, double ycen, double fwhm, int xbin,
              int ybin, int ndiv, bool comp_dfwhm) {

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

  double *dsky_ptr = static_cast<double *>(dsky_info.ptr);
  double *dheight_ptr = static_cast<double *>(dheight_info.ptr);
  double *dxcen_ptr = static_cast<double *>(dxcen_info.ptr);
  double *dycen_ptr = static_cast<double *>(dycen_info.ptr);

  py::buffer_info dfwhm_info;
  double *dfwhm_ptr = nullptr;

  if (comp_dfwhm) {
    dfwhm_info = dfwhm.request();
    dfwhm_ptr = static_cast<double *>(dfwhm_info.ptr);
  }

  // Calculate Gaussian profile parameter
  double alpha = 4.0 * std::log(2.0) / (fwhm * fwhm);
  double two_alpha_height = 2.0 * alpha * height;
  double dfwhm_coeff = two_alpha_height / fwhm;

  size_t n_pixels = x_info.shape[0] * x_info.shape[1];

  // Initialize dsky to ones (derivative of sky is always 1)
  std::fill_n(dsky_ptr, n_pixels, 1.0);

  if (ndiv > 0) {
    // With sub-pixellation
    double inv_nadd = 1.0 / static_cast<double>(xbin * ybin * ndiv * ndiv);
    double inv_ndiv = 1.0 / static_cast<double>(ndiv);

    // Mean offset within sub-pixels
    double soff = (ndiv - 1.0) / (2.0 * ndiv);

    // Loop over all pixels
    for (size_t pixel_idx = 0; pixel_idx < n_pixels; ++pixel_idx) {
      double x_val = x_ptr[pixel_idx];
      double y_val = y_ptr[pixel_idx];
      double dheight = 0.0;
      double dxcen = 0.0;
      double dycen = 0.0;
      double dfwhm = 0.0;

      // Loop over sub-pixels
      for (int iy = 0; iy < ybin; ++iy) {
        double yoff = iy - (ybin - 1) / 2.0 - soff;
        for (int ix = 0; ix < xbin; ++ix) {
          double xoff = ix - (xbin - 1) / 2.0 - soff;
          for (int isy = 0; isy < ndiv; ++isy) {
            double ysoff = yoff + isy * inv_ndiv;
            for (int isx = 0; isx < ndiv; ++isx) {
              double xsoff = xoff + isx * inv_ndiv;
              double dx = x_val + xsoff - xcen;
              double dy = y_val + ysoff - ycen;
              double rsq = dx * dx + dy * dy;

              // Gaussian value
              double dh = std::exp(-alpha * rsq);
              dheight += dh;
              dxcen += two_alpha_height * dh * dx;
              dycen += two_alpha_height * dh * dy;

              if (comp_dfwhm) {
                dfwhm += dfwhm_coeff * dh * rsq;
              }
            }
          }
        }
      }

      dheight_ptr[pixel_idx] = dheight * inv_nadd;
      dxcen_ptr[pixel_idx] = dxcen * inv_nadd;
      dycen_ptr[pixel_idx] = dycen * inv_nadd;

      if (comp_dfwhm) {
        dfwhm_ptr[pixel_idx] = dfwhm * inv_nadd;
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

  m.def("gaussian", &gaussian_cpp, "C++ implementation of Gaussian profile",
        py::arg("x"), py::arg("y"), py::arg("sky"), py::arg("height"),
        py::arg("xcen"), py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"),
        py::arg("ybin"), py::arg("ndiv"));

  m.def("dgaussian", &dgaussian_cpp,
        "C++ implementation of Gaussian profile derivatives", py::arg("x"),
        py::arg("y"), py::arg("sky"), py::arg("height"), py::arg("xcen"),
        py::arg("ycen"), py::arg("fwhm"), py::arg("xbin"), py::arg("ybin"),
        py::arg("ndiv"), py::arg("comp_dfwhm"));
}
