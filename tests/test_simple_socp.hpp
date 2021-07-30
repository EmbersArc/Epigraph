#include <array>

using namespace cvx;

// This example solves a simple random second order cone problem
// based on https://www.cvxpy.org/examples/basic/socp.html

TEST_CASE("Simple SOCP")
{
    // Set up problem data.

    // ! This test's feasibility depends on the rng.
    srand(123);

    // number of second order cone constraints
    const size_t m = 3;
    // number of variables
    const size_t n = 10;
    // dimension of equality constraints
    const size_t p = 5;
    // dimension of second order cone constraints
    const size_t n_i = 5;

    std::array<Eigen::Matrix<double, n_i, n>, m> A;
    std::array<Eigen::Matrix<double, n_i, 1>, m> b;
    std::array<Eigen::Matrix<double, n, 1>, m> c;
    std::array<double, m> d;

    Eigen::Matrix<double, n, 1> x0;
    x0.setRandom();
    Eigen::Matrix<double, n, 1> f;
    f.setRandom();

    for (size_t i = 0; i < m; i++)
    {
        A[i].setRandom();
        b[i].setRandom();
        c[i].setRandom();
        d[i] = (A[i] * x0).norm() - c[i].dot(x0);
    }

    Eigen::Matrix<double, p, n> F;
    F.setRandom();
    Eigen::Matrix<double, p, 1> g = F * x0;

    // Formulate SOCP.
    auto t0 = std::chrono::high_resolution_clock::now();

    // Create the SOCP instance.
    OptimizationProblem socp;

    // Add variables. Those can be scalars, vectors or matrices.
    VectorX x = socp.addVariable("x", n);

    // Add constraints.
    // SOCP
    for (size_t i = 0; i < m; i++)
    {
        socp.addConstraint(lessThan((par(A[i]) * x + par(b[i])).norm(),
                                    par(c[i]).dot(x) + par(d[i])));
    }
    // Equality
    socp.addConstraint(equalTo(par(F) * x, par(g)));

    // Here we use dynamic parameter. This allows changing it without reformulating the problem.
    socp.addCostTerm(dynpar(f).transpose() * x);

    // Print the problem for inspection.
    std::cout << socp << "\n";

    // Create the solver instance.
    ecos::ECOSSolver solver(socp);

    // Solve the problem and show solver output.
    t0 = std::chrono::high_resolution_clock::now();
    const bool success = solver.solve(false);
    REQUIRE(solver.isFeasible(1e-8));

    if (not success)
    {
        // This should not happen in this example.
        throw std::runtime_error("Solver returned a critical error.");
    }

    std::cout << "Solver message: " << solver.getResultString() << "\n";

    // Get Solution.
    Eigen::Matrix<double, n, 1> x_eval = eval(x);
    // Print the first solution.
    std::cout << "First solution:\n"
              << x_eval << "\n\n";
    // Test the constraints
    for (size_t i = 0; i < m; i++)
    {
        REQUIRE((A[i] * x_eval + b[i]).norm() <= Approx(c[i].dot(x_eval) + d[i]).margin(1e-6));
    }
    REQUIRE((F * x_eval - g).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-8));

    // Change the problem parameters and solve again.
    f.setRandom();
    solver.solve(false);
    REQUIRE(solver.isFeasible(1e-8));
    x_eval = eval(x);

    // Test the constraints again
    for (size_t i = 0; i < m; i++)
    {
        REQUIRE((A[i] * x_eval + b[i]).norm() <= Approx(c[i].dot(x_eval) + d[i]).margin(1e-6));
    }
    REQUIRE((F * x_eval - g).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-8));

    std::cout << "Solution after changing the cost function:\n"
              << x_eval << "\n\n";
}