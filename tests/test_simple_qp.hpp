// A very basic QP example

TEST_CASE("Simple QP 1")
{
    // Set up problem data.
    Eigen::MatrixXd P(2, 2);
    Eigen::VectorXd q(2);
    Eigen::MatrixXd A(3, 2);
    Eigen::VectorXd l(3);
    Eigen::VectorXd u(3);

    P << 2, 0.5, 0.5, 1;
    q << 1, 1;
    l << 1, 0, 0;
    u << 1, 0.7, 0.7;
    A << 1, 1, 1, 0, 0, 1;

    // Formulate QP.
    OptimizationProblem qp;

    VectorX x = var("x", 2);

    qp.addConstraint(box(par(l), par(A) * x, par(u)));

    qp.addCostTerm(x.transpose() * par(P) * x);
    qp.addCostTerm(par(q).transpose() * x);

    // Create and initialize the solver instance.
    osqp::OSQPSolver solver(qp);
    solver.setAlpha(1.0);

    // Solve the problem and show solver output.
    solver.solve(true);

    fmt::print("{}", qp);
    fmt::print("{}", solver);

    Eigen::VectorXd x_val = eval(x);
    Eigen::Vector2d x_sol(0.3, 0.7);

    const double cost = x_val.transpose() * P * x_val + q.dot(x_val);

    fmt::print("Solution:\n {}\n", x_val);
    fmt::print("Cost:    \n {}\n", cost);

    REQUIRE(x_val.isApprox(x_sol, 1e-5));

    // ! The order in which the variables are added to the problem may vary, so this test fails randomly.
    // auto test_stream = std::ostringstream();
    // test_stream << solver;
    // REQUIRE(test_stream.str() == "Quadratic problem\nMinimize 0.5x'Px + q'x\nSubject to l <= Ax <= u\nWith:\n\nP:\n4 1\n1 2\n\nq:\n1\n1\n\nA:\n1 1\n1 0\n0 1\n\nl:\n1\n0\n0\n\nu:\n  1\n0.7\n0.7");
}

TEST_CASE("Simple QP 2")
{
    OptimizationProblem qp;

    VectorX x = var("x", 3);
    qp.addConstraint(equalTo(x.sum(), 1.));
    qp.addConstraint(box(-1., x, 1.));
    qp.addCostTerm((2. + x(1)) * x(1) + (1. + x(0)) * x(0) + (1. + x(0)) * x(1) + x(2) * (2. + x(2)) + x(2) * x(2));

    osqp::OSQPSolver solver(qp);

    fmt::print("{}\n", qp);
    fmt::print("{}\n", solver);

    solver.solve(true);

    Eigen::Vector3d x_eval = eval(x);

    REQUIRE((x_eval - Eigen::Vector3d(1., -1. / 3., 1. / 3.)).cwiseAbs().maxCoeff() <= 1e-3);
}

TEST_CASE("Non-convex QP")
{
    OptimizationProblem qp;

    VectorX x = var("x", 3);

    qp.addConstraint(equalTo(x.sum(), 1.));
    qp.addConstraint(box(-1., x, 1.));

    Eigen::Matrix3d M;
    M.setZero();
    M.diagonal() << -3, -2, -1;

    qp.addCostTerm(x.transpose() * par(M) * x);

    fmt::print("{}\n", qp);

    REQUIRE_THROWS(osqp::OSQPSolver(qp));
}
