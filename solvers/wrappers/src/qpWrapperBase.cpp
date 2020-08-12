#include "qpWrapperBase.hpp"

#include <Eigen/Cholesky>

namespace cvx::internal
{

    QPWrapperBase::QPWrapperBase(OptimizationProblem &problem)
    {
        std::vector<Eigen::Triplet<Parameter>> A_coeffs, P_coeffs;
        std::vector<Parameter> l_coeffs, u_coeffs;

        // Build equality constraint parameters
        for (internal::EqualityConstraint &constraint : problem.equality_constraints)
        {
            constraint.affine.cleanUp();
            if (constraint.affine.isConstant())
            {
                continue;
            }

            for (Term &term : constraint.affine.terms)
            {
                addVariable(term.variable);
                A_coeffs.emplace_back(u_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }
            l_coeffs.push_back(Parameter(-1.) * constraint.affine.constant);
            u_coeffs.push_back(Parameter(-1.) * constraint.affine.constant);
        }

        // Build positive constraint parameters
        for (internal::PositiveConstraint &constraint : problem.positive_constraints)
        {
            constraint.affine.cleanUp();
            if (constraint.affine.isConstant())
            {
                continue;
            }

            for (Term &term : constraint.affine.terms)
            {
                addVariable(term.variable);
                A_coeffs.emplace_back(u_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }
            l_coeffs.push_back(Parameter(-1.) * constraint.affine.constant);
            u_coeffs.push_back(Parameter(std::numeric_limits<double>::max()));
        }

        // Build box constraint parameters
        for (internal::BoxConstraint &constraint : problem.box_constraints)
        {
            if (constraint.lower.isConstant() and constraint.upper.isConstant())
            {
                // lower <= middle <= upper
                constraint.middle.cleanUp();
                if (constraint.middle.isConstant())
                {
                    continue;
                }

                for (Term &term : constraint.middle.terms)
                {
                    addVariable(term.variable);
                    A_coeffs.emplace_back(u_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                l_coeffs.push_back(constraint.lower.constant - constraint.middle.constant);
                u_coeffs.push_back(constraint.upper.constant - constraint.middle.constant);
            }
            else
            {
                // c_lower - c_middle <= middle - lower <= inf
                Affine middle_m_lower = constraint.middle - constraint.lower;
                middle_m_lower.cleanUp();

                if (middle_m_lower.isFirstOrder())
                {
                    for (Term &term : middle_m_lower.terms)
                    {
                        addVariable(term.variable);
                        A_coeffs.emplace_back(u_coeffs.size(),
                                              term.variable.getProblemIndex(),
                                              term.parameter);
                    }
                    l_coeffs.push_back(constraint.lower.constant - constraint.middle.constant);
                    u_coeffs.push_back(Parameter(std::numeric_limits<double>::max()));
                }

                // c_middle - c_upper <= upper - middle <= inf
                Affine upper_m_middle = constraint.upper - constraint.middle;
                upper_m_middle.cleanUp();

                if (upper_m_middle.isFirstOrder())
                {
                    for (Term &term : upper_m_middle.terms)
                    {
                        addVariable(term.variable);
                        A_coeffs.emplace_back(u_coeffs.size(),
                                              term.variable.getProblemIndex(),
                                              term.parameter);
                    }
                    l_coeffs.push_back(constraint.middle.constant - constraint.upper.constant);
                    u_coeffs.push_back(Parameter(std::numeric_limits<double>::max()));
                }
            }
        }

        // Build cost function
        if (problem.costFunction.getOrder() == 0 or problem.costFunction.isNorm())
        {
            throw std::runtime_error("QP cost functions must be linear or quadratic.");
        }

        // Linear part
        for (Term &term : problem.costFunction.affine.terms)
        {
            q_params(term.variable.getProblemIndex()) += term.parameter;
        }

        // Quadratic part
        for (Product &product : problem.costFunction.products)
        {
            for (Term &term1 : product.firstTerm().terms)
            {
                for (Term &term2 : product.secondTerm().terms)
                {
                    addVariable(term1.variable);
                    addVariable(term2.variable);

                    const std::pair<size_t, size_t> sorted = std::minmax(term1.variable.getProblemIndex(),
                                                                         term2.variable.getProblemIndex());

                    Parameter param = term1.parameter * term2.parameter;

                    // Explicitly double diagonal elements
                    if (sorted.first == sorted.second)
                    {
                        param *= Parameter(2.);
                    }

                    P_coeffs.emplace_back(sorted.first,
                                          sorted.second,
                                          param);
                }
            }

            // The linear parts from the multiplication
            if (not product.firstTerm().constant.isZero())
            {
                for (Term &term : product.secondTerm().terms)
                {
                    addVariable(term.variable);
                    q_params(term.variable.getProblemIndex()) += product.firstTerm().constant * term.parameter;
                }
            }
            if (not product.secondTerm().constant.isZero())
            {
                for (Term &term : product.firstTerm().terms)
                {
                    addVariable(term.variable);
                    q_params(term.variable.getProblemIndex()) += product.secondTerm().constant * term.parameter;
                }
            }
        }

        assert(l_coeffs.size() == u_coeffs.size());

        // Fill matrices and vectors
        A_params.resize(l_coeffs.size(), getNumVariables());
        P_params.resize(getNumVariables(), getNumVariables());

        A_params.setFromTriplets(A_coeffs.begin(), A_coeffs.end());
        P_params.setFromTriplets(P_coeffs.begin(), P_coeffs.end());

        l_params = Eigen::Map<VectorXp>(l_coeffs.data(), l_coeffs.size());
        u_params = Eigen::Map<VectorXp>(u_coeffs.data(), u_coeffs.size());

        solution.resize(getNumVariables());
    }

    size_t QPWrapperBase::getNumInequalityConstraints() const
    {
        return A_params.rows();
    }

    void QPWrapperBase::addVariable(Variable &variable)
    {
        if (not variable.isLinkedToSolver())
        {
            variable.linkToSolver(&solution, getNumVariables());
            variables.push_back(variable);
            q_params.conservativeResize(getNumVariables());
        }
    }

    std::ostream &operator<<(std::ostream &os, const QPWrapperBase &wrapper)
    {
        Eigen::MatrixXd A = eval(wrapper.A_params);
        Eigen::MatrixXd P = eval(wrapper.P_params);
        P += P.triangularView<Eigen::StrictlyUpper>().transpose();
        Eigen::VectorXd q = eval(wrapper.q_params);
        Eigen::VectorXd l = eval(wrapper.l_params);
        Eigen::VectorXd u = eval(wrapper.u_params);

        os << "Quadratic problem\n";
        os << "Minimize 0.5x'Px + q'x\n";
        os << "Subject to l <= Ax <= u\n";
        os << "With:\n\n";

        os << "P:\n"
           << P << "\n\n";
        os << "q:\n"
           << q << "\n\n";
        os << "A:\n"
           << A << "\n\n";
        os << "l:\n"
           << l << "\n\n";
        os << "u:\n"
           << u;

        return os;
    }

    bool QPWrapperBase::isConvex() const
    {
        if (P_params.nonZeros() == 0)
        {
            return true;
        }
        else
        {
            Eigen::LLT<Eigen::MatrixXd> llt(eval(P_params));
            return llt.info() != Eigen::NumericalIssue;
        }
    }

} // namespace cvx
