// This example solves the portfolio optimization problem in QP form

using namespace cvx;

TEST_CASE("Portfolio SOCP")
{
    const size_t n = 5; // Assets
    const size_t m = 2; // Factors

    // Set up problem data.
    double gamma = 0.5;          // risk aversion parameter
    Eigen::VectorXd mu(n);       // vector of expected returns
    Eigen::MatrixXd F(n, m);     // factor-loading matrix
    Eigen::VectorXd D(n);        // diagonal of idiosyncratic risk
    Eigen::MatrixXd Sigma(n, n); // asset return covariance

    // mu.setRandom();
    // F.setRandom();
    // D.setRandom();
    // mu = mu.cwiseAbs();
    // F = F.cwiseAbs();
    // D = D.cwiseAbs();
    // Sigma = F * F.transpose();
    // Sigma.diagonal() += D;

    mu << 0.680375, 0.211234, 0.566198, 0.59688, 0.823295;
    F << 0.604897, 0.0452059, 0.329554, 0.257742, 0.536459, 0.270431, 0.444451, 0.0268018, 0.10794, 0.904459;
    Sigma << 1.20033, 0.210998, 0.336728, 0.270059, 0.106179, 0.210998, 0.44646, 0.246494, 0.153379, 0.268689, 0.336728, 0.246494, 0.795515, 0.245678, 0.302499, 0.270059, 0.153379, 0.245678, 0.91505, 0.0722151, 0.106179, 0.268689, 0.302499, 0.0722151, 1.04364;
    D << 0.83239, 0.271423, 0.434594, 0.716795, 0.213938;

    // Formulate SOCP.
    OptimizationProblem socp;

    VectorX x = socp.addVariable("x", n);
    Scalar t = socp.addVariable("t");
    Scalar s = socp.addVariable("s");
    Scalar u = socp.addVariable("u");
    Scalar v = socp.addVariable("v");

    socp.addConstraint(greaterThan(x, 0.));
    socp.addConstraint(equalTo(x.sum(), 1.));
    socp.addConstraint(lessThan((par(D.cwiseSqrt()).cwiseProduct(x)).norm(), u));
    socp.addConstraint(lessThan((par(F).transpose() * x).norm(), v));
    VectorX lhs1(2);
    lhs1 << 1. - t, 2. * u;
    socp.addConstraint(lessThan(lhs1.norm(), 1. + t));
    VectorX lhs2(2);
    lhs2 << 1. - s, 2. * v;
    socp.addConstraint(lessThan(lhs2.norm(), 1. + s));

    socp.addCostTerm(gamma * (t + s) - dynpar(mu).dot(x));

    // Create the solver instance.
    cvx::ecos::ECOSSolver solver(socp);

    solver.solve(false);
    {
        const Eigen::VectorXd x_eval = eval(x);
        Eigen::VectorXd x_sol(n);
        x_sol << 0.24424712, 0., 0.01413456, 0.25067381, 0.4909445;
        REQUIRE((x_eval - x_sol).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-4));
        REQUIRE(x_eval.minCoeff() >= 0.);
        REQUIRE(x_eval.sum() == Approx(1.));
    }

    // Update data
    // mu.setRandom();
    // mu = mu.cwiseAbs();
    mu << 0.967399, 0.514226, 0.725537, 0.608354, 0.686642;

    // Solve again
    solver.solve(false);
    {
        const Eigen::VectorXd x_eval = eval(x);
        Eigen::VectorXd x_sol(n);
        x_sol << 4.38579051e-01, 3.04375987e-23, 2.00025310e-01, 1.17002001e-01, 2.44393639e-01;
        REQUIRE((x_eval - x_sol).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-4));
        REQUIRE(x_eval.minCoeff() >= 0.);
        REQUIRE(x_eval.sum() == Approx(1.));
        REQUIRE(solver.isFeasible(1e-8));
    }
}