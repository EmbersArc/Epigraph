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

        auto construct = [&op]() { eicos::EiCOSSolver solver(op); };

        REQUIRE_THROWS(construct());
    }
}
