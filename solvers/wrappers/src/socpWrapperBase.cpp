#include "socpWrapperBase.hpp"

namespace cvx
{

    SOCPWrapperBase::SOCPWrapperBase(OptimizationProblem &problem)
    {
        std::vector<Eigen::Triplet<Parameter>> A_coeffs, G_coeffs;
        std::vector<Parameter> b_coeffs, h_coeffs;
        std::vector<int> cone_dimensions;

        // Build equality constraint parameters (b - A * x == 0)
        for (EqualityConstraint &c : problem.equality_constraints)
        {
            for (Term &term : c.affine.terms)
            {
                addVariable(term.variable);
                A_coeffs.emplace_back(b_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }

            b_coeffs.push_back(c.affine.constant);
        }

        // Build positive constraint parameters
        for (PositiveConstraint &constraint : problem.positive_constraints)
        {
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
        for (BoxConstraint &constraint : problem.box_constraints)
        {
            // lower <= middle <= upper

            // 0 <= middle - lower
            for (Term &term : constraint.middle.terms)
            {
                addVariable(term.variable);
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }
            for (Term &term : constraint.lower.terms)
            {
                addVariable(term.variable);
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      -term.parameter);
            }
            h_coeffs.push_back(constraint.middle.constant - constraint.lower.constant);

            // 0 <= upper - middle
            for (Term &term : constraint.upper.terms)
            {
                addVariable(term.variable);
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      term.parameter);
            }
            for (Term &term : constraint.middle.terms)
            {
                addVariable(term.variable);
                G_coeffs.emplace_back(h_coeffs.size(),
                                      term.variable.getProblemIndex(),
                                      -term.parameter);
            }
            h_coeffs.push_back(constraint.upper.constant - constraint.middle.constant);
        }

        // Build second order cone constraint parameters
        for (SecondOrderConeConstraint &constraint : problem.second_order_cone_constraints)
        {
            // Affine part
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
        if (problem.costFunction.getOrder() != 1)
        {
            throw std::runtime_error("SOCP cost functions must be linear.");
        }

        c_params.resize(n_variables);
        c_params.setZero();

        for (Term &term : problem.costFunction.affine.terms)
        {
            if (term.variable.isLinkedToProblem())
                c_params(term.variable.getProblemIndex()) += term.parameter;
        }

        // Fill matrices and vectors
        A_params.resize(b_coeffs.size(), n_variables);
        G_params.resize(h_coeffs.size(), n_variables);

        A_params.setFromTriplets(A_coeffs.begin(), A_coeffs.end());
        G_params.setFromTriplets(G_coeffs.begin(), G_coeffs.end());
        b_params = Eigen::Map<VectorXp>(b_coeffs.data(), b_coeffs.size());
        h_params = Eigen::Map<VectorXp>(h_coeffs.data(), h_coeffs.size());
        soc_dims = Eigen::Map<Eigen::VectorXi>(cone_dimensions.data(), cone_dimensions.size());

        solution.resize(n_variables);
    }

    size_t SOCPWrapperBase::getNumVariables() const
    {
        return n_variables;
    }

    void SOCPWrapperBase::addVariable(Variable &variable)
    {
        if (not variable.isLinkedToProblem())
        {
            variable.linkToProblem(&solution, n_variables);
            n_variables++;
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
