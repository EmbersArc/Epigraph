#include "wrappers/socpWrapperBase.hpp"

namespace cvx::internal
{

    SOCPWrapperBase::SOCPWrapperBase(OptimizationProblem &problem)
    {
        std::vector<Eigen::Triplet<Parameter>> A_coeffs, G_coeffs;
        std::vector<Parameter> b_coeffs, h_coeffs;
        std::vector<int> cone_dimensions;

        // Build equality constraint parameters (b - A * x == 0)
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
                A_coeffs.emplace_back(b_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }

            b_coeffs.push_back(constraint.affine.constant);
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
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }

            h_coeffs.push_back(constraint.affine.constant);
        }

        // Build box constraint parameters
        for (internal::BoxConstraint &constraint : problem.box_constraints)
        {
            // lower <= middle <= upper

            // 0 <= middle - lower
            Affine middle_m_lower = constraint.middle - constraint.lower;
            middle_m_lower.cleanUp();

            if (middle_m_lower.isFirstOrder())
            {
                for (Term &term : middle_m_lower.terms)
                {
                    addVariable(term.variable);
                    G_coeffs.emplace_back(h_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                h_coeffs.push_back(middle_m_lower.constant);
            }

            // 0 <= upper - middle
            Affine upper_m_middle = constraint.upper - constraint.middle;
            upper_m_middle.cleanUp();

            if (upper_m_middle.isFirstOrder())
            {
                for (Term &term : upper_m_middle.terms)
                {
                    addVariable(term.variable);
                    G_coeffs.emplace_back(h_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                h_coeffs.push_back(upper_m_middle.constant);
            }
        }

        // Build second order cone constraint parameters
        for (internal::SecondOrderConeConstraint &constraint : problem.second_order_cone_constraints)
        {
            // Affine part
            constraint.affine.cleanUp();

            for (Term &term : constraint.affine.terms)
            {
                addVariable(term.variable);
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }
            h_coeffs.push_back(constraint.affine.constant);

            // Norm part
            for (Affine &affine : constraint.norm)
            {
                affine.cleanUp();

                if (affine.isZero())
                {
                    continue;
                }

                for (Term &term : affine.terms)
                {
                    addVariable(term.variable);
                    G_coeffs.emplace_back(h_coeffs.size(),
                                          term.variable.getProblemIndex(),
                                          term.parameter);
                }
                h_coeffs.push_back(affine.constant);
            }

            cone_dimensions.push_back(constraint.norm.size() + 1);
        }

        // Build cost function
        problem.costFunction.affine.cleanUp();
        if (problem.costFunction.getOrder() != 1)
        {
            throw std::runtime_error("SOCP cost functions must be linear.");
        }

        for (Term &term : problem.costFunction.affine.terms)
        {
            addVariable(term.variable);
            c_params(term.variable.getProblemIndex()) += term.parameter;
        }

        // Fill matrices and vectors
        A_params.resize(b_coeffs.size(), getNumVariables());
        G_params.resize(h_coeffs.size(), getNumVariables());

        A_params.setFromTriplets(A_coeffs.begin(), A_coeffs.end());
        G_params.setFromTriplets(G_coeffs.begin(), G_coeffs.end());
        b_params = Eigen::Map<VectorXp>(b_coeffs.data(), b_coeffs.size());
        h_params = Eigen::Map<VectorXp>(h_coeffs.data(), h_coeffs.size());
        soc_dims = Eigen::Map<Eigen::VectorXi>(cone_dimensions.data(), cone_dimensions.size());

        solution->resize(getNumVariables());
    }

    void SOCPWrapperBase::addVariable(Variable &variable)
    {
        const bool was_linked = variable.linkToSolver(solution, getNumVariables());
        if (was_linked)
        {
            variables.push_back(variable);
            c_params.conservativeResize(getNumVariables());
        }
    }

    size_t SOCPWrapperBase::getNumEqualityConstraints() const
    {
        return A_params.rows();
    }

    size_t SOCPWrapperBase::getNumInequalityConstraints() const
    {
        return G_params.rows();
    }

    size_t SOCPWrapperBase::getNumPositiveConstraints() const
    {
        return G_params.rows() - soc_dims.sum();
    }

    size_t SOCPWrapperBase::getNumCones() const
    {
        return soc_dims.size();
    }

    std::ostream &operator<<(std::ostream &os, const SOCPWrapperBase &wrapper)
    {
        Eigen::VectorXd c = eval(wrapper.c_params);
        Eigen::MatrixXd G = -eval(wrapper.G_params);
        Eigen::VectorXd h = eval(wrapper.h_params);
        Eigen::MatrixXd A = -eval(wrapper.A_params);
        Eigen::VectorXd b = eval(wrapper.b_params);

        os << "Second order cone problem\n";
        os << "Minimize c'x\n";
        os << "Subject to Gx <=_K h\n";
        os << "           Ax == b\n";
        os << "With:\n\n";

        os << "c:\n"
           << c << "\n\n";
        os << "G:\n"
           << G << "\n\n";
        os << "h:\n"
           << h << "\n\n";
        os << "A:\n"
           << A << "\n\n";
        os << "b:\n"
           << b << "\n\n";

        return os;
    }

} // namespace cvx
