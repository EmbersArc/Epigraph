using namespace cvx;

TEST_CASE("Simple Problems")
{
    {
        OptimizationProblem op;
        VectorX x = var("x", 2);

        REQUIRE(eval(x(0)) == 0.);
        REQUIRE(eval(x(1)) == 0.);

        op.addConstraint(greaterThan(x, 1.));
        op.addConstraint(lessThan(x.norm(), 5.));
        op.addCostTerm(-x.sum());

        eicos::EiCOSSolver solver(op);

        solver.solve(true);

        Eigen::VectorXd x_eval = eval(x);
        const double optval_eval = solver.getInfo().pcost;

        const double x_sol = std::sqrt(25. / 2.);
        const double optval_sol = -2 * x_sol;

        REQUIRE(x_sol == Approx(x_eval(0)).margin(1e-7));
        REQUIRE(x_sol == Approx(x_eval(1)).margin(1e-7));
        REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
    }
    {
        OptimizationProblem op;
        VectorX x = var("x", 2);

        op.addConstraint(greaterThan(x, 1.));
        op.addConstraint(lessThan(sqrt(x(0) * x(0) + x(1) * x(1) + par(2.)), 5.));
        op.addCostTerm(-x.sum());

        eicos::EiCOSSolver solver(op);

        solver.solve(true);

        Eigen::VectorXd x_eval = eval(x);
        const double optval_eval = solver.getInfo().pcost;

        const double x_sol = std::sqrt((5. * 5. - 2.) / 2.);
        const double optval_sol = -2 * x_sol;

        REQUIRE(x_sol == Approx(x_eval(0)).margin(1e-7));
        REQUIRE(x_sol == Approx(x_eval(1)).margin(1e-7));
        REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
    }
    { // Testing the edge case with (p1 * x1) * (p2 * x1)
        OptimizationProblem op;
        VectorX x = var("x", 2);

        op.addConstraint(greaterThan(x, 1.));
        Scalar lhs = sqrt((par(2.) * x(0)) * (par(3.) * x(0)) + x(1) * x(1));
        op.addConstraint(lessThan(lhs, 5.));
        op.addCostTerm(-x.sum());

        eicos::EiCOSSolver solver(op);

        solver.solve(true);

        Eigen::VectorXd x_eval = eval(x);
        const double optval_eval = solver.getInfo().pcost;

        Eigen::VectorXd x_sol(2);
        x_sol << 1., 4.35889894;
        const double optval_sol = -x_sol.sum();

        REQUIRE(x_sol(0) == Approx(x_eval(0)).margin(1e-7));
        REQUIRE(x_sol(1) == Approx(x_eval(1)).margin(1e-7));
        REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
    }
    { // Invalid problem
        OptimizationProblem op;
        VectorX x = var("x", 2);

        op.addCostTerm(x(0) * x(1));

        REQUIRE_THROWS(eicos::EiCOSSolver(op));
    }
}

TEST_CASE("Least Squares")
{
    Eigen::MatrixXd A(20, 15);
    Eigen::VectorXd b(20);
    Eigen::VectorXd x_sol(15);

    A << 1.62434536, -0., -0., -0., 0.,
        -2.3015387, 1.74481176, -0., 0.3190391, -0.,
        1.46210794, -2.06014071, -0.3224172, -0., 1.13376944,
        -1.09989127, -0.17242821, -0., 0., 0.,
        -0., 0., 0.90159072, 0., 0.90085595,
        -0.68372786, -0., -0.93576943, -0., 0.53035547,
        -0.69166075, -0.39675353, -0., -0., -0.67124613,
        -0., -0., 0., 0., 0.,
        -0., -0.88762896, -0., 1.6924546, 0.,
        -0.63699565, 0., 0., 0.12015895, 0.,
        0.30017032, -0.35224985, -0., -0.34934272, -0.,
        0., 0.83898341, 0., 0., 0.,
        -0., 1.25286816, 0., -0.29809284, 0.,
        -0.07557171, 1.13162939, 0., 2.18557541, -1.39649634,
        -1.44411381, -0., 0.16003707, 0., 0.31563495,
        -2.02220122, -0.30620401, 0., 0., 0.76201118,
        -0.22232814, -0.20075807, 0., 0., 0.19829972,
        0., -0., 0., 0., 0.,
        0., 0.18515642, -0., -0.63873041, 0.,
        0.07734007, -0., 0.04359686, -0., 0.69803203,
        -0., 1.2245077, 0.40349164, 0.59357852, -0.,
        0.16938243, 0., -0., -0.26621851, 0.03261455,
        -0., 0.31515939, 0.84616065, -0., 0.,
        -0., -0., -1.61577235, 1.12141771, 0.,
        -0., -0.77516162, 0., 0., -1.85798186,
        1.23616403, 0., 0.3380117, -1.19926803, 0.,
        -0.1809203, -0.60392063, -1.23005814, 0., 0.79280687,
        -0., 0.52057634, -1.14434139, 0.80186103, 0.0465673,
        -0.18656977, -0., 0., 0.75041164, 0.,
        0.13770121, 0., 0.61838026, 0., 0.68255141,
        -0., -2.43483776, 1.0388246, 0., 0.44136444,
        -0., -0., -0.11905419, 0., -0.,
        -0., -0., 0.24879916, -0.29664115, 0.,
        -0., 0.98633519, 0.2135339, 2.19069973, -0.,
        -0.64691669, 0., 0., -0.24863478, 0.,
        -0., 0., -0., 0.68006984, -0.,
        -1.27255876, 0., 0., 1.29322588, -0.11044703,
        -0.61736206, 0.5627611, 0.24073709, 0.28066508, -0.,
        0., 0.36949272, 0., 1.1110567, 0.,
        -0., 0.60231928, 0., 0., 0.,
        -0., 0.82400562, -0., 0., -0.,
        -1.76068856, -0., -0.89055558, -0., 1.9560789,
        -0.3264995, -1.34267579, 1.11438298, -0., -1.23685338,
        0., 0., -0.43495668, 1.40754, 0.12910158,
        1.6169496, 0., 1.55880554, 0.1094027, -1.2197444,
        2.44936865, -0.54577417, -0., -0.7003985, -0.20339445,
        0.24266944, 0., 0., 1.79215821, -0.,
        -0., -1.18231813, -0., -0., 0.82502982,
        -0., -0., -0., -1.39662042, -0.,
        0., 0., -0.44317193, 1.81053491, -0.,
        -0., -0., -2.793085, 0., 0.,
        -1.04458938, 0., 0.585662, 0., -0.6069984,
        0.10622272, -0., 0.79502609, -0.37443832, 0.,
        1.20205486, 0.28474811, 0.26246745, 0., -0.7332716,
        0., 1.54335911, 0.75880566, 0.88490881, -0.87728152,
        -0.86778722, -1.44087602, 0., -0., 1.39984394,
        -0.78191168, -0., 0., 0., 0.0607502,
        0., 0.01652757, 0., -1.11647002, 0.,
        -0., -0., 0.49233656, -0.68067814, -0.08450803,
        -0.29736188, 0.417302, 0.78477065, -0., 0.58591043;

    b << -1.07296428, 0.49515861, -0.9520621, -0.51814555, -1.4614036,
        -0.51634791, 0.3511169, -0.06877046, -1.34776494, 1.47073986,
        0.33722094, 1.00806543, 0.78522692, -0.66486777, -1.94504696,
        -0.91542437, 1.22515585, -1.05354607, 0.81604368, -0.61240697;

    x_sol << 0.39859857,
        -0.76549,
        -1.58716462,
        0.53706734,
        -0.07398277,
        -1.29907978,
        -0.66583103,
        0.07945784,
        0.06242338,
        0.17445139,
        -0.60960874,
        0.9141227,
        -0.5453288,
        -0.24262929,
        -0.67958804;

    {
        OptimizationProblem qp;

        VectorX x = var("x", 15);

        qp.addCostTerm((par(A) * x - par(b)).squaredNorm());

        fmt::print("{}\n", qp);

        osqp::OSQPSolver solver(qp);
        solver.solve(true);
        REQUIRE((eval(x) - x_sol).cwiseAbs().maxCoeff() < 1e-5);
    }
    {
        OptimizationProblem qp;

        VectorX x = var("x", 15);
        Scalar s = var("s");

        // Since norms are equivalent
        qp.addConstraint(lessThan((par(A) * x - par(b)).norm(), s));
        qp.addCostTerm(s);

        fmt::print("{}\n", qp);

        eicos::EiCOSSolver solver(qp);
        solver.solve(true);
        REQUIRE((eval(x) - x_sol).cwiseAbs().maxCoeff() < 1e-5);
    }
}