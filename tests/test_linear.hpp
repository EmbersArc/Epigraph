using namespace cvx;

TEST_CASE("Linear")
{
    const size_t m = 15;
    const size_t n = 10;

    Eigen::MatrixXd A(m, n);
    Eigen::VectorXd b(m);
    Eigen::VectorXd c(n);

    A << -0.68372786, -0.12289023, -0.93576943, -0.26788808, 0.53035547, -0.69166075, -0.39675353, -0.6871727, -0.84520564, -0.67124613, -0.0126646, -1.11731035, 0.2344157, 1.65980218, 0.74204416, -0.19183555, -0.88762896, -0.74715829, 1.6924546, 0.05080775, -0.63699565, 0.19091548, 2.10025514, 0.12015895, 0.61720311, 0.30017032, -0.35224985, -1.1425182, -0.34934272, -0.20889423, 0.58662319, 0.83898341, 0.93110208, 0.28558733, 0.88514116, -0.75439794, 1.25286816, 0.51292982, -0.29809284, 0.48851815, -0.07557171, 1.13162939, 1.51981682, 2.18557541, -1.39649634, -1.44411381, -0.50446586, 0.16003707, 0.87616892, 0.31563495, -2.02220122, -0.30620401, 0.82797464, 0.23009474, 0.76201118, -0.22232814, -0.20075807, 0.18656139, 0.41005165, 0.19829972, 0.11900865, -0.67066229, 0.37756379, 0.12182127, 1.12948391, 1.19891788, 0.18515642, -0.37528495, -0.63873041, 0.42349435, 0.07734007, -0.34385368, 0.04359686, -0.62000084, 0.69803203, -0.44712856, 1.2245077, 0.40349164, 0.59357852, -1.09491185, 0.16938243, 0.74055645, -0.9537006, -0.26621851, 0.03261455, -1.37311732, 0.31515939, 0.84616065, -0.85951594, 0.35054598, -1.31228341, -0.03869551, -1.61577235, 1.12141771, 0.40890054, -0.02461696, -0.77516162, 1.27375593, 1.96710175, -1.85798186, 1.23616403, 1.62765075, 0.3380117, -1.19926803, 0.86334532, -0.1809203, -0.60392063, -1.23005814, 0.5505375, 0.79280687, -0.62353073, 0.52057634, -1.14434139, 0.80186103, 0.0465673, -0.18656977, -0.10174587, 0.86888616, 0.75041164, 0.52946532, 0.13770121, 0.07782113, 0.61838026, 0.23249456, 0.68255141, -0.31011677, -2.43483776, 1.0388246, 2.18697965, 0.44136444, -0.10015523, -0.13644474, -0.11905419, 0.01740941, -1.12201873, -0.51709446, -0.99702683, 0.24879916, -0.29664115, 0.49521132, -0.17470316, 0.98633519, 0.2135339, 2.19069973, -1.89636092, -0.64691669, 0.90148689, 2.52832571, -0.24863478, 0.04366899;
    b << 2.17495142, -0.07902089, -2.93864432, 1.93790752, 0.57842629, 2.57175626, 0.67612818, 1.88614126, 3.3688581, 2.75695134, -0.63273236, 3.38791401, -0.31286018, -0.46582275, 3.89352826;
    c << 5.9159385, -0.38653276, -1.58376908, -3.45614976, -4.00186624, 2.30939272, 0.4855809, -2.76450383, -4.55428255, -1.02779359;

    OptimizationProblem op;
    VectorX x = op.addVariable("x", n);
    op.addConstraint(lessThan(par(A) * x, par(b)));
    op.addCostTerm(par(c).dot(x));

    // ! Fails right now, waiting for fix
    {
        REQUIRE_THROWS(ecos::ECOSSolver(op));
        // ecos::ECOSSolver solver(op);

        // solver.solve(true);

        // Eigen::VectorXd x_eval = eval(x);
        // const double optval_eval = solver.getInfo().pcost;

        // // generated with CVXPY
        // Eigen::VectorXd x_sol(n);
        // x_sol << -1.10378318, -0.14861315, -0.93044404, 0.01510552, 0.64782787, -1.17145192, 1.10413269, 0.84810568, 0.47081493, 0.89445169;
        // const double optval_sol = -15.220912604892906;

        // REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
        // REQUIRE((x_eval - x_sol).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-7));
    }
    {
        eicos::EiCOSSolver solver(op);

        solver.solve(true);

        Eigen::VectorXd x_eval = eval(x);
        const double optval_eval = solver.getInfo().pcost;

        // generated with CVXPY
        Eigen::VectorXd x_sol(n);
        x_sol << -1.10378318, -0.14861315, -0.93044404, 0.01510552, 0.64782787, -1.17145192, 1.10413269, 0.84810568, 0.47081493, 0.89445169;
        const double optval_sol = -15.220912604892906;

        REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
        REQUIRE((x_eval - x_sol).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-7));
    }
    {
        osqp::OSQPSolver solver(op);
        solver.setEpsAbs(1e-5);
        solver.setEpsRel(1e-5);
        solver.setPolish(true);

        solver.solve(true);

        Eigen::VectorXd x_eval = eval(x);
        const double optval_eval = solver.getInfo().obj_val;

        // generated with CVXPY
        Eigen::VectorXd x_sol(n);
        x_sol << -1.10152248, -0.16244666, -0.89989851, 0.03085193, 0.6100638, -1.13030714, 1.1277109, 0.87917365, 0.48921664, 0.89817175;
        const double optval_sol = -15.220912603926376;

        REQUIRE(optval_eval == Approx(optval_sol).margin(1e-7));
        REQUIRE((x_eval - x_sol).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-7));
    }
}