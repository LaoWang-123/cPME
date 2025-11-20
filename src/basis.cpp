// src/basis.cpp
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

// -----------------------------------------------------------------------------
// Helper: compute cpq normalization constant
// -----------------------------------------------------------------------------
inline double cpq_val(int p, int q) {
  if (p == 0 && q == 0) {
    return 1.0;
  } else if (p == 0 || q == 0) {
    return std::sqrt(2.0);
  } else {
    return 2.0;
  }
}

//
// ============================================================================
// (1) Scalar Fourier basis ψ_{pq}(u, v)
// ============================================================================
//

//' Scalar Fourier basis \eqn{\psi_{pq}(u,v)}.
//'
//' Compute the scalar Fourier basis function
//' \deqn{\psi_{pq}(u,v) = c_{pq}\cos(2\pi p u)\cos(2\pi q v).}
//'
//' @param u Numeric value in \[0,1\].
//' @param v Numeric value in \[0,1\].
//' @param p Integer Fourier index in the u-direction.
//' @param q Integer Fourier index in the v-direction.
//'
//' @return A numeric scalar.
//' @export
// [[Rcpp::export]]
double scalar_basis_cpp(double u, double v, int p, int q) {
 double cpq = cpq_val(p, q);
 double A   = 2.0 * M_PI * static_cast<double>(p);
 double B   = 2.0 * M_PI * static_cast<double>(q);

 return cpq * std::cos(A * u) * std::cos(B * v);
}

//
// ============================================================================
// (2) Gradient of Fourier basis ∇ψ_{pq}(u, v)
// ============================================================================
//

//' Gradient of Fourier basis \eqn{\nabla\psi_{pq}(u,v)}.
//'
//' Computes the gradient:
//' \deqn{
//'   \nabla\psi_{pq}(u,v) = \left(
//'     -c_{pq}(2\pi p)\sin(2\pi p u)\cos(2\pi q v),\;
//'     -c_{pq}(2\pi q)\cos(2\pi p u)\sin(2\pi q v)
//'   \right).
//' }
//'
//' @param u Numeric input.
//' @param v Numeric input.
//' @param p Integer Fourier index.
//' @param q Integer Fourier index.
//'
//' @return A numeric vector of length 2.
//' @export
// [[Rcpp::export]]
NumericVector grad_basis_cpp(double u, double v, int p, int q) {
 double cpq = cpq_val(p, q);
 double A   = 2.0 * M_PI * static_cast<double>(p);
 double B   = 2.0 * M_PI * static_cast<double>(q);

 NumericVector out(2);

 out[0] = -cpq * A * std::sin(A * u) * std::cos(B * v);
 out[1] = -cpq * B * std::cos(A * u) * std::sin(B * v);

 return out;
}

//
// ============================================================================
// (3) Rotated gradient *∇ψ_{pq}(u, v)
// ============================================================================
//

//' Rotated gradient \eqn{*\nabla\psi_{pq}(u,v)}.
//'
//' Computes the 90-degree rotated gradient:
//' \deqn{
//'   *\nabla\psi_{pq}(u,v) = \left(
//'     c_{pq}(2\pi q)\cos(2\pi p u)\sin(2\pi q v),\;
//'    -c_{pq}(2\pi p)\sin(2\pi p u)\cos(2\pi q v)
//'   \right).
//' }
//'
//' @param u Numeric input.
//' @param v Numeric input.
//' @param p Integer Fourier index.
//' @param q Integer Fourier index.
//'
//' @return A numeric vector of length 2.
//' @export
// [[Rcpp::export]]
NumericVector rot_basis_cpp(double u, double v, int p, int q) {
 double cpq = cpq_val(p, q);
 double A   = 2.0 * M_PI * static_cast<double>(p);
 double B   = 2.0 * M_PI * static_cast<double>(q);

 NumericVector out(2);

 out[0] =  cpq * B * std::cos(A * u) * std::sin(B * v);
 out[1] = -cpq * A * std::sin(A * u) * std::cos(B * v);

 return out;
}

//
// ============================================================================
// (4) Jacobian of gradient D(∇ψ_{pq})
// ============================================================================
//

//' Jacobian of gradient field \eqn{D(\nabla\psi_{pq})(u,v)}.
//'
//' Computes the 2×2 Jacobian matrix of the gradient field.
//'
//' @param u Numeric input.
//' @param v Numeric input.
//' @param p Integer Fourier index.
//' @param q Integer Fourier index.
//'
//' @return A 2×2 numeric matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix jacobian_grad_cpp(double u, double v, int p, int q) {
 double cpq = cpq_val(p, q);
 double A   = 2.0 * M_PI * static_cast<double>(p);
 double B   = 2.0 * M_PI * static_cast<double>(q);

 NumericMatrix J(2, 2);

 double cosAu = std::cos(A * u);
 double cosBv = std::cos(B * v);
 double sinAu = std::sin(A * u);
 double sinBv = std::sin(B * v);

 J(0, 0) = -A * A * cpq * cosAu * cosBv;
 J(0, 1) =  A * B * cpq * sinAu * sinBv;
 J(1, 0) =  A * B * cpq * sinAu * sinBv;
 J(1, 1) = -B * B * cpq * cosAu * cosBv;

 return J;
}

//
// ============================================================================
// (5) Jacobian of rotated gradient D(*∇ψ_{pq})
// ============================================================================
//

//' Jacobian of rotated gradient field \eqn{D(*\nabla\psi_{pq})(u,v)}.
//'
//' Computes the 2×2 Jacobian matrix of the rotated gradient.
//'
//' @param u Numeric input.
//' @param v Numeric input.
//' @param p Integer Fourier index.
//' @param q Integer Fourier index.
//'
//' @return A 2×2 numeric matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix jacobian_rot_cpp(double u, double v, int p, int q) {
 double cpq = cpq_val(p, q);
 double A   = 2.0 * M_PI * static_cast<double>(p);
 double B   = 2.0 * M_PI * static_cast<double>(q);

 NumericMatrix J(2, 2);

 double cosAu = std::cos(A * u);
 double cosBv = std::cos(B * v);
 double sinAu = std::sin(A * u);
 double sinBv = std::sin(B * v);

 J(0, 0) = -A * B * cpq * sinAu * sinBv;
 J(0, 1) =  B * B * cpq * cosAu * cosBv;
 J(1, 0) = -A * A * cpq * cosAu * cosBv;
 J(1, 1) =  A * B * cpq * sinAu * sinBv;

 return J;
}
