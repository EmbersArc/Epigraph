TEST_CASE("getVariableValue")
{
    OptimizationProblem qp;

    double s_eval;
    Eigen::VectorXd v_eval;
    Eigen::MatrixXd m_eval;

    Scalar s = qp.addVariable("s");
    VectorX v = qp.addVariable("v", 3);
    MatrixX m = qp.addVariable("m", 3, 3);

    REQUIRE_THROWS(qp.getVariableValue("imaginary_s", s_eval));
    REQUIRE_THROWS(qp.getVariableValue("imaginary_v", v_eval));
    REQUIRE_THROWS(qp.getVariableValue("imaginary_m", m_eval));

    qp.getVariableValue("s", s_eval);
    qp.getVariableValue("v", v_eval);
    qp.getVariableValue("m", m_eval);

    REQUIRE(s_eval == 0.);
    REQUIRE(v_eval == Eigen::Vector3d::Zero());
    REQUIRE(m_eval == Eigen::Matrix3d::Zero());
}

TEST_CASE("getProblemIndex")
{
    Variable variable("x");
    REQUIRE_THROWS(variable.getProblemIndex());
}

TEST_CASE("getVariable")
{
    OptimizationProblem op;
    Scalar scalar = op.addVariable("scalar");
    VectorX vector = op.addVariable("vector", 5);
    MatrixX matrix = op.addVariable("matrix", 5, 5);

    REQUIRE_THROWS(op.addVariable("scalar"));
    REQUIRE_THROWS(op.addVariable("vector", 5));
    REQUIRE_THROWS(op.addVariable("matrix", 5, 5));

    Scalar scalar_returned;
    VectorX vector_returned;
    MatrixX matrix_returned;

    op.getVariable("scalar", scalar_returned);
    op.getVariable("vector", vector_returned);
    op.getVariable("matrix", matrix_returned);

    REQUIRE(scalar == scalar_returned);
    REQUIRE(vector == vector_returned);
    REQUIRE(matrix == matrix_returned);

    REQUIRE_THROWS(op.getVariable("imaginary_scalar", scalar_returned));
    REQUIRE_THROWS(op.getVariable("imaginary_vector", vector_returned));
    REQUIRE_THROWS(op.getVariable("imaginary_matrix", matrix_returned));
}
