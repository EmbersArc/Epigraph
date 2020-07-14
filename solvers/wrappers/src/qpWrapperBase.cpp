#include "qpWrapperBase.hpp"

#include <Eigen/Cholesky>

namespace cvx
{

    QPWrapperBase::QPWrapperBase(OptimizationProblem &problem)
    {
        std::vector<Eigen::Triplet<Parameter>> A_coeffs, P_coeffs;
        std::vector<Parameter> l_coeffs, u_coeffs;

        // Build equality constraint parameters
        for (EqualityConstraint &constraint : problem.equality_constraints)
        {
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
        for (PositiveConstraint &constraint : problem.positive_constraints)
        {
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
        for (BoxConstraint &constraint : problem.box_constraints)
        {
            if (constraint.lower.isConstant() and constraint.upper.isConstant())
            {
                // lower <= middle <= upper
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
                for (Term &term : constraint.middle.terms)
                {
                    addVariable(term.variable);
                    A_coeffs.emplace_back(u_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                for (Term &term : constraint.lower.terms)
                {
                    addVariable(term.variable);
                    A_coeffs.emplace_back(u_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          -term.parameter);
                }
                l_coeffs.push_back(constraint.lower.constant - constraint.middle.constant);
                u_coeffs.push_back(Parameter(std::numeric_limits<double>::max()));

                // c_middle - c_upper <= upper - middle <= inf
                for (Term &term : constraint.upper.terms)
                {
                    addVariable(term.variable);
                    A_coeffs.emplace_back(u_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                for (Term &term : constraint.middle.terms)
                {
                    addVariable(term.variable);
                    A_coeffs.emplace_back(u_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          -term.parameter);
                }
                l_coeffs.push_back(constraint.middle.constant - constraint.upper.constant);
                u_coeffs.push_back(Parameter(std::numeric_limits<double>::max()));
            }
        }

        // Build cost function
        if (problem.costFunction.isNorm())
        {
            throw std::runtime_error("QP cost functions must be linear or quadratic.");
        }

        // Linear part
        q_params.resize(n_variables);
        q_params.setZero();
        for (Term &term : problem.costFunction.affine.terms)
        {
            if (term.variable.isLinkedToProblem())
                q_params(term.variable.getProblemIndex()) += term.parameter;
        }

        // Quadratic part
        for (Product &product : problem.costFunction.products)
        {
            for (Term &term1 : product.firstTerm().terms)
            {
                for (Term &term2 : product.secondTerm().terms)
                {
                    // Only relevant when both variables are constrained. Otherwise, one can be zero.
                    if (term1.variable.isLinkedToProblem() and term2.variable.isLinkedToProblem())
                    {
                        const std::pair<size_t, size_t> sorted = std::minmax(term1.variable.getProblemIndex(),
                                                                             term2.variable.getProblemIndex());

                        Parameter param = term1.parameter * term2.parameter;
                        if (sorted.first == sorted.second)
                            param *= Parameter(2.);

                        P_coeffs.emplace_back(sorted.first,
                                              sorted.second,
                                              param);
                    }
                }
            }

            // The linear parts from the multiplication
            if (not product.firstTerm().constant.isZero())
            {
                for (Term &term : product.secondTerm().terms)
                {
                    if (term.variable.isLinkedToProblem())
                    {
                        q_params(term.variable.getProblemIndex()) += product.firstTerm().constant * term.parameter;
                    }
                }
            }
            if (not product.secondTerm().constant.isZero())
            {
                for (Term &term : product.firstTerm().terms)
                {
                    if (term.variable.isLinkedToProblem())
                    {
                        q_params(term.variable.getProblemIndex()) += product.secondTerm().constant * term.parameter;
                    }
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

        solution.resize(n_variables);
    }

    size_t QPWrapperBase::getNumVariables() const
    {
        return n_variables;
    }

    size_t QPWrapperBase::getNumInequalityConstraints() const
    {
        return A_params.rows();
    }

    void QPWrapperBase::addVariable(Variable &variable)
    {
        if (not variable.isLinkedToProblem())
        {
            variable.linkToProblem(&solution, n_variables);
            n_variables++;
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
        Eigen::LLT<Eigen::MatrixXd> llt(eval(P_params));
        return not(llt.info() == Eigen::NumericalIssue);
    }

} // namespace cvx
