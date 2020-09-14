using namespace cvx;

TEST_CASE("Constraint")
{
    {
        OptimizationProblem op;

        VectorX x = op.addVariable("x", 2);

        auto test_stream = std::ostringstream();
        test_stream << equalTo(x.sum(), 1.);
        REQUIRE(test_stream.str() == "x[0] + x[1] + -1 == 0");

        test_stream = std::ostringstream();
        test_stream << lessThan(x(0), 1.);
        REQUIRE(test_stream.str() == "0 <= -1 * x[0] + 1");

        test_stream = std::ostringstream();
        test_stream << greaterThan(1., x(0));
        REQUIRE(test_stream.str() == "0 <= -1 * x[0] + 1");

        test_stream = std::ostringstream();
        test_stream << lessThan(x.norm(), 1.);
        REQUIRE(test_stream.str() == "((x[0])^2 + (x[1])^2)^(1/2) <= 1");

        test_stream = std::ostringstream();
        test_stream << box(-1., x.sum(), 1.);
        REQUIRE(test_stream.str() == "-1 <= x[0] + x[1] <= 1");
    }

    { // Test box constraint with multiple variables
        OptimizationProblem op;

        VectorX x = op.addVariable("x", 6);
        Scalar lhs = x(0) + x(1) + 1.;
        Scalar mid = x(2) + x(3);
        Scalar rhs = x(4) + x(5) + 1.;

        op.addConstraint(greaterThan(x, 0.));
        op.addConstraint(box(lhs, mid, rhs));
        op.addCostTerm(par(Eigen::VectorXd::Random(6).cwiseAbs()).dot(x));

        {
            osqp::OSQPSolver solver(op);

            solver.solve(false);

            REQUIRE(eval(lhs) <= eval(mid) + 1e-3);
            REQUIRE(eval(mid) <= eval(rhs) + 1e-3);
            REQUIRE(op.getOptimalValue() == Approx(solver.getInfo().obj_val));
        }

        {
            ecos::ECOSSolver solver(op);

            solver.solve(false);

            REQUIRE(eval(lhs) <= eval(mid) + 1e-3);
            REQUIRE(eval(mid) <= eval(rhs) + 1e-3);
            REQUIRE(op.getOptimalValue() == Approx(solver.getInfo().pcost));
        }
    }
    { // Invalid constraints
        OptimizationProblem op;
        VectorX x = op.addVariable("x", 3);

        REQUIRE_THROWS(lessThan(x(0), x.norm()));
        REQUIRE_THROWS(lessThan(x.squaredNorm(), x(0)));
        REQUIRE_THROWS(box(x(0), x.squaredNorm(), x(1)));
        REQUIRE_THROWS(equalTo(x, x.squaredNorm()));
    }
}