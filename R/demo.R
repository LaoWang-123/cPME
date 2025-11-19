### Run Demo

grad_f_bowl <- function(u, v, R = 1.0, c = 1.5) {
    # ensure vectorized operation
    du <- 2 * R
    dv <- 2 * R
    dz_du <- 4 * c * (2 * u - 1)
    dz_dv <- 4 * c * (2 * v - 1)
    
    # 返回一个 3×2 矩阵
    matrix(c(
      du, 0,
      0,  dv,
      dz_du, dz_dv
    ), nrow = 3, byrow = TRUE)
}

f2_grad_fn <- grad_f_bowl

basis_set <- build_basis_set(3,3)

Ugrid <- expand.grid(
  u = seq(0, 1, by = 0.01),
  v = seq(0, 1, by = 0.01)
)

u=0.4
v=0.6
