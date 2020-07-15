using namespace cvx;

#include <iostream>
#include <string>

TEST_CASE("Scalar")
{
    {
        VectorX x = var("x", 2);
        REQUIRE_NOTHROW(x.dot(x) + x.dot(x));
        REQUIRE_THROWS(x.norm() + x.norm());
        REQUIRE_THROWS(x.norm() + x.dot(x));
        REQUIRE_THROWS(x.dot(x) + x.norm());
        REQUIRE_THROWS(x.dot(x) - x.dot(x));
        REQUIRE_THROWS(x.dot(x) * x.dot(x));
        REQUIRE_THROWS(x(0) / x.sum());
        REQUIRE_THROWS(x.squaredNorm() / par(2.));
        REQUIRE_THROWS(sqrt(x.sum()));

        auto test_stream = std::ostringstream();
        test_stream << x.sum();
        REQUIRE(test_stream.str() == "x[0] + x[1]");

        test_stream = std::ostringstream();
        test_stream << par(2.) * x.sum();
        REQUIRE(test_stream.str() == "2 * x[0] + 2 * x[1]");

        test_stream = std::ostringstream();
        Scalar test_scalar = x.norm() + par(1.);
        test_stream << test_scalar;
        REQUIRE(test_stream.str() == "((x[0])^2 + (x[1])^2)^(1/2) + 1");

        test_stream = std::ostringstream();
        test_scalar = x(0) * x(0) + x(0) * x(1) + x(0) + par(1.);
        test_stream << test_scalar;
        REQUIRE(test_stream.str() == "(x[0])^2 + (x[0]) * (x[1]) + x[0] + 1");

        // Check the equality operator
        REQUIRE(test_scalar == x(0) * x(0) + x(0) * x(1) + x(0) + par(1.));
    }

    {
        Scalar scalar = var("s");
        VectorX vector = var("v", 2);
        MatrixX matrix = var("m", 2, 2);

        auto test_stream = std::ostringstream();
        test_stream << scalar;
        REQUIRE(test_stream.str() == "s");

        test_stream = std::ostringstream();
        test_stream << vector.sum();
        REQUIRE(test_stream.str() == "v[0] + v[1]");

        test_stream = std::ostringstream();
        test_stream << matrix.col(0).sum();
        REQUIRE(test_stream.str() == "m[0, 0] + m[1, 0]");
    }

    {
        VectorX x = var("x", 3);

        OptimizationProblem op;

        op.addCostTerm(x.sum());
        op.addConstraint(greaterThan(x, par(Eigen::Vector3d(1., 2., 3.))));

        osqp::OSQPSolver solver(op);

        solver.setAlpha(1.);

        solver.solve();

        Eigen::Vector3d x_sol(1., 2., 3.);

        REQUIRE(eval(x.norm()) == Approx(x_sol.norm()));
        REQUIRE(eval(x.norm() + x.sum()) == Approx(x_sol.norm() + x_sol.sum()));
        REQUIRE(eval(x(0) * x(1)) == Approx(x_sol(0) * x_sol(1)));
        REQUIRE(x(0) * x(1) == x(1) * x(0));
        REQUIRE(eval(x(0) / par(2.)) == Approx(x_sol(0) / 2.));
    }

    { // Dynamic parameters
        double d1 = 1.;
        double d2 = 2.;

        Scalar p1 = dynpar(d1);
        Scalar p2 = dynpar(d2);

        REQUIRE(eval(p1 * p2) == d1 * d2);

        d1 = 2.;
        d2 = 3.;

        REQUIRE(eval(p1 * p2) == d1 * d2);
    }
}