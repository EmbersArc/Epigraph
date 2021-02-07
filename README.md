<p align="center">
    <img src="img/logo.png">
</p>

<p align="center">
<a href="https://shields.io/" alt="">
    <img src="https://img.shields.io/badge/Version-0.4.1-orange.svg" /></a>
<a href="https://shields.io/" alt="">
    <img src="https://img.shields.io/badge/Status-Beta-orange.svg" /></a>
<a href="https://shields.io/" alt="">
    <img src="https://img.shields.io/badge/Maintained%3F-yes-green.svg" /></a>
<a href="https://github.com/EmbersArc/Epigraph/blob/master/LICENSE" alt="">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" /></a>
<a href="https://travis-ci.org/github/EmbersArc/Epigraph" alt="">
    <img src="https://api.travis-ci.org/EmbersArc/Epigraph.svg?branch=master" /></a>
<a href="https://codecov.io/gh/EmbersArc/Epigraph" alt="">
    <img src="https://codecov.io/gh/EmbersArc/Epigraph/branch/master/graph/badge.svg" /></a>
</p>

Epigraph is a modern C++ interface to formulate and solve linear, quadratic and second order cone problems. It makes use of Eigen types and operator overloading for straightforward problem formulation.

## Features
* Flexible and intuitive way to formulate LPs, QPs and SOCPs
* Dynamic parameters that can be changed without re-formulating the problem
* Automatically clean up the problem and remove unused variables
* Print the problem formulation and solver data for inspection

## Dependencies

* [Eigen](http://eigen.tuxfamily.org)

## Supported Solvers

The solvers are included as submodules for convenience. Note that some solvers have more restrictive licenses which automatically override the Epigraph license when activated. Pass the listed argument to cmake during configuration to enable the solvers.

### QP 
* [OSQP](https://github.com/oxfordcontrol/osqp) `-DENABLE_OSQP=TRUE`. Apache-2.0 License.

### SOCP
* [ECOS](https://github.com/embotech/ecos) `-DENABLE_ECOS=TRUE`. GPLv3 License.

## Usage

### CMake
To use Epigraph with a cmake project, simply enable the desired solvers, include the subdirectory and link the library.
```
set(ENABLE_OSQP TRUE)
set(ENABLE_ECOS TRUE)
add_subdirectory(Epigraph)
target_link_libraries(my_library epigraph)
```

### Documentation

While the example below is likely enough to get you started, more explanation is provided below. The full reference can be found [here](https://embersarc.github.io/Epigraph/).

### Example

```cpp
#include "epigraph.hpp"

#include <iostream>

// This example solves the portfolio optimization problem in QP form

using namespace cvx;

int main()
{
    size_t n = 5; // Assets
    size_t m = 2; // Factors

    // Set up problem data.
    double gamma = 0.5;          // risk aversion parameter
    Eigen::VectorXd mu(n);       // vector of expected returns
    Eigen::MatrixXd F(n, m);     // factor-loading matrix
    Eigen::VectorXd D(n);        // diagonal of idiosyncratic risk
    Eigen::MatrixXd Sigma(n, n); // asset return covariance

    mu.setRandom();
    F.setRandom();
    D.setRandom();

    mu = mu.cwiseAbs();
    F = F.cwiseAbs();
    D = D.cwiseAbs();
    Sigma = F * F.transpose();
    Sigma.diagonal() += D;

    // Formulate QP.
    OptimizationProblem qp;

    // Declare variables with...
    // addVariable(name) for scalars,
    // addVariable(name, rows) for vectors and
    // addVariable(name, rows, cols) for matrices.
    VectorX x = qp.addVariable("x", n);

    // Available constraint types are equalTo(), lessThan(), greaterThan() and box()
    qp.addConstraint(greaterThan(x, 0.));
    qp.addConstraint(equalTo(x.sum(), 1.));

    // Make mu dynamic in the cost function so we can change it later
    qp.addCostTerm(x.transpose() * par(gamma * Sigma) * x - dynpar(mu).dot(x));

    // Print the problem formulation for inspection
    std::cout << qp << "\n";

    // Create and initialize the solver instance.
    osqp::OSQPSolver solver(qp);

    // Print the canonical problem formulation for inspection
    std::cout << solver << "\n";

    // Solve problem and show solver output
    const bool verbose = true;
    solver.solve(verbose);

    std::cout << "Solver message:  " << solver.getResultString() << "\n";
    std::cout << "Solver exitcode: " << solver.getExitCode() << "\n";

    // Call eval() to get the variable values
    std::cout << "Solution:\n" << eval(x) << "\n";

    // Update data
    mu.setRandom();
    mu = mu.cwiseAbs();

    // Solve again
    // OSQP will warm start automatically
    solver.solve(verbose);

    std::cout << "Solution after changing the cost function:\n" << eval(x) << "\n";
}
```
See the [tests](tests) for more examples, including the same problem in SOCP form.

### Variables
Variables are created by directly adding them to a problem:
```cpp
    OptimizationProblem qp;
    cvx::Scalar scalar_var = qp.addVariable("x");
    cvx::VectorX vector_var = qp.addVariable("x", n);
    cvx::MatrixX matrix_var = qp.addVariable("x", n, m);
```
They contain the solution values after the problem has been solved successfully. Those values can be retrieved with the `eval` function, or by casting the value to a `double`.
```cpp
    double scalar_sol = cvx::eval(scalar_var);
    Eigen::VectorXd vector_sol = cvx::eval(vector_var);
    Eigen::MatrixXd matrix_sol = cvx::eval(matrix_var);
```

### Parameters
There are three kinds of parameters:
#### Constant
A value that can't be changed after instantiating the solver. Use the `par` function to turn scalars or Eigen types into constant parameters.
#### Dynamic
A value that *can* be changed after instantiating the solver. Use the `dynpar` function to turn scalars or Eigen types into dynamic parameters. Internally, this stores a pointer to the original data and will fetch the data each time the problem is solved. Important: Do not move or let this data go out of scope before the solver instance.
#### Operation
This parameter type is created when using the operations `+`, `-`, `*` or `/` with dynamic parameters. This records the operations and will later execute them again to build the new problem based on the changed dynamic parameters. Using said operations with constant parameters will again yield constant parameters and not result in any additional computations.

### Problem Formulation
The following terms may be passed to the constraint functions:

| Function | Allowed expressions |
| --- | --- |
| `equalTo()`|`Affine == Affine` |
| `lessThan()`| `Affine <= Affine` or `Norm2 + Affine <= Affine` (SOCP) |
| `greaterThan()`| `Affine >= Affine` or `Affine >= Norm2 + Affine` (SOCP) |
| `box()`| `Affine <= Affine <= Affine` |
| `addCostTerm()`| `Affine` (SOCP) or `QuadForm + Affine` (QP) |

With the following expressions:

| Expression | Form |
| --- | --- |
| `Affine` | `p1 * x1 + p2 * x2 + ... + c` |
| `Norm2` | `(Affine1^2  + Affine2^2 + ...)^(1/2)` |
| `QuadForm` | ``x' * P * x`` where `P` is Hermitian |
