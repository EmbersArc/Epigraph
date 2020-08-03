
using namespace cvx;

TEST_CASE("MPC QP")
{
    // This test solves an MPC problem with quadratic programming.

    size_t T = 7; // time horizon

    Eigen::MatrixXd A(2, 2);
    A << 2, -1, 1, 0.2;
    Eigen::MatrixXd B(2, 1);
    B << 1, 0;
    Eigen::VectorXd x0(2);
    x0 << 3, 1;

    OptimizationProblem qp;

    // Create variables
    MatrixX x = qp.addVariable("x", 2, T + 1);
    MatrixX u = qp.addVariable("u", 1, T);

    // Dynamics
    for (size_t t = 0; t < T; t++)
    {
        qp.addConstraint(equalTo(x.col(t + 1), par(A) * x.col(t) + par(B) * u.col(t)));
    }

    // State and control limits
    qp.addConstraint(box(-5., x, 5.));

    // Equivalent to
    // qp.addConstraint(box(-2., u, 2.));
    qp.addConstraint(greaterThan(u, -2.));
    qp.addConstraint(lessThan(u, 2.));

    // Boundary constraints
    qp.addConstraint(equalTo(x.col(0), par(x0)));
    qp.addConstraint(equalTo(x.col(T), 0.));

    // Cost function
    qp.addCostTerm(x.squaredNorm() + u.squaredNorm());

    // Print the problem for inspection
    fmt::print("{}", qp);

    // Create and initialize the solver instance.
    osqp::OSQPSolver solver(qp);

    solver.setAlpha(1.);

    // Solve the problem and get the solution
    solver.solve(true);
    fmt::print("Solver result: {} ({})\n", solver.getResultString(), solver.getExitCode());
    Eigen::MatrixXd x_sol = eval(x);
    Eigen::MatrixXd u_sol = eval(u);

    fmt::print("X:\n {}\n\n", x_sol);
    fmt::print("U:\n {}\n\n", u_sol);

    // Check solution
    for (size_t t = 0; t < T; t++)
    {
        const Eigen::VectorXd x_solution = x_sol.col(t + 1);
        const Eigen::VectorXd x_propagate = A * x_sol.col(t) + B * u_sol.col(t);
        const double max_error = (x_propagate - x_solution).cwiseAbs().maxCoeff();
        REQUIRE(max_error == Approx(0.).margin(1e-5));
    }
    REQUIRE(x_sol.maxCoeff() <= Approx(5.).margin(1e-3));
    REQUIRE(x_sol.minCoeff() >= Approx(-5.).margin(1e-3));
    REQUIRE(u_sol.maxCoeff() <= Approx(2.).margin(1e-3));
    REQUIRE(u_sol.minCoeff() >= Approx(-2.).margin(1e-3));
}