# C-Cox: High-Performance Cox Proportional Hazards Regression

C-Cox is a lightweight, standalone C implementation of the Cox Proportional Hazards (PH) model. It is designed for high performance survival analysis on large datasets where memory efficiency and computational speed are critical.

The model estimates the effects of covariates on the hazard function, defined as:

$$\lambda(t | X) = \lambda_0(t) \exp\left(\sum_{i=1}^{p} \beta_i X_i\right)$$

---

## Why C-Cox?

* **High Performance:** Optimized Newton-Raphson solver for rapid convergence.
* **Efficient Risk Set Calculation:** Reduced computational complexity for partial likelihood gradients.
* **Minimal Dependencies:** Written in standard C with no requirement for heavy external math libraries.
* **Tie Handling:** Supports standard methods including Breslow and Efron approximations.
* **Portability:** Easily integrated into existing C/C++ projects or wrapped for Python/R/Julia via FFI.

---

## Some specs:

* **Language:** C99 or higher.
* **Optimization:** Iterative Partial Likelihood Maximization.
* **Memory Profile:** Fixed size allocations suited for high concurrency environments.
* **Input Format:** Row major or Column major matrix support for survival data (time, status, covariates).

---

## Installation

To compile the library, use any standard C compiler (e.g., GCC or Clang):

```bash
sudo apt update && sudo apt install -y libgsl-dev
make
./cox_model config.dat
