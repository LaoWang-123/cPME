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

// ============================================================================
// (1) Scalar Fourier basis ψ_{pq}(u, v)
// ============================================================================

/**
 * @title Scalar Fourier basis ψ_{pq}(u, v)
 * @name scalar_basis_cpp
 * @description
 * Compute the scalar Fourier basis function
 * \f[
 *   \psi_{pq}(u, v) = c_{pq} \cos(2\pi p u) \cos(2\pi q v).
 * \f]
 *
 * @param u A numeric value in [0,1].
 * @param v A numeric value in [0,1].
 * @param p Integer Fourier index along the u-direction.
 * @param q Integer Fourier index along the v-direction.
 *
 * @return A numeric scalar, the evaluated basis value.
 * @export
 */
// [[Rcpp::export]]
double scalar_basis_cpp(double u, double v, int p, int q) {
  double cpq = cpq_val(p, q);
  double A   = 2.0 * M_PI * static_cast<double>(p);
  double B   = 2.0 * M_PI * static_cast<double>(q);

  return cpq * std::cos(A * u) * std::cos(B * v);
}

// ============================================================================
// (2) Gradient of scalar basis ∇ψ_{pq}(u, v)
// ============================================================================

/**
 * @title Gradient of Fourier basis ∇ψ_{pq}(u, v)
 * @name grad_basis_cpp
 * @description
 * Compute the gradient vector field of the Fourier basis function:
 * \f[
 *   \nabla \psi_{pq}(u,v)
 *   = \left(
 *       -c_{pq} (2\pi p) \sin(2\pi p u)\cos(2\pi q v),\;
 *       -c_{pq} (2\pi q) \cos(2\pi p u)\sin(2\pi q v)
 *     \right).
 * \f]
 *
 * @param u A numeric value.
 * @param v A numeric value.
 * @param p Integer Fourier index p.
 * @param q Integer Fourier index q.
 *
 * @return A numeric vector of length 2 representing (du-component, dv-component).
 * @export
 */
// [[Rcpp::export]]
NumericVector grad_basis_cpp(double u, double v, int p, int q) {
  double cpq = cpq_val(p, q);
  double A   = 2.0 * M_PI * static_cast<double>(p);
  double B   = 2.0 * M_PI * static_cast<double>(q);

  NumericVector out(2);

  out[0] = -cpq * A * std::sin(A * u) * std::cos(B * v); // ∂ψ/∂u
  out[1] = -cpq * B * std::cos(A * u) * std::sin(B * v); // ∂ψ/∂v

  return out;
}

// ============================================================================
// (3) Rotated gradient *∇ψ_{pq}(u, v)
// ============================================================================

/**
 * @title Rotated gradient of Fourier basis *∇ψ_{pq}(u, v)
 * @name rot_basis_cpp
 * @description
 * Compute the 90-degree rotated gradient field
 * \f[
 *   *\nabla \psi_{pq}(u,v)
 *   = \left(
 *       c_{pq} (2\pi q)\cos(2\pi p u)\sin(2\pi q v),\;
 *      -c_{pq} (2\pi p)\sin(2\pi p u)\cos(2\pi q v)
 *     \right).
 * \f]
 *
 * @param u numeric input.
 * @param v numeric input.
 * @param p Fourier index p.
 * @param q Fourier index q.
 *
 * @return A numeric vector of length 2.
 * @export
 */
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

// ============================================================================
// (4) Jacobian of gradient basis D(∇ψ_{pq})
// ============================================================================

/**
 * @title Jacobian of gradient D(∇ψ_{pq})(u,v)
 * @name jacobian_grad_cpp
 * @description
 * Compute the 2×2 Jacobian matrix of the gradient field of ψ_{pq}.
 *
 * Rows correspond to derivatives w.r.t (u, v), and columns correspond to
 * gradient components.
 *
 * @param u numeric input.
 * @param v numeric input.
 * @param p Fourier index p.
 * @param q Fourier index q.
 *
 * @return A 2×2 numeric matrix.
 * @export
 */
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

// ============================================================================
// (5) Jacobian of rotated gradient D(*∇ψ_{pq})
// ============================================================================

/**
 * @title Jacobian of rotated gradient D(*∇ψ_{pq})(u,v)
 * @name jacobian_rot_cpp
 * @description
 * Compute the 2×2 Jacobian matrix of the rotated gradient field *∇ψ_{pq}.
 *
 * @param u numeric input.
 * @param v numeric input.
 * @param p Fourier index p.
 * @param q Fourier index q.
 *
 * @return A 2×2 numeric matrix.
 * @export
 */
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
