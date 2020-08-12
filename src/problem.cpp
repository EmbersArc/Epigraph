#include "problem.hpp"

namespace cvx
{
    using namespace internal;

    Scalar OptimizationProblem::addVariable(const std::string &name)
    {
        if (scalar_variables.find(name) != scalar_variables.end())
        {
            const std::string error_message = "Could not add scalar variable '" + name + "' since it already exists.";
            throw std::runtime_error(error_message);
        }

        Scalar variable = Variable(name);

        scalar_variables.emplace(name, variable);

        return variable;
    }

    VectorX OptimizationProblem::addVariable(const std::string &name,
                                             size_t rows)
    {
        if (vector_variables.find(name) != vector_variables.end())
        {
            const std::string error_message = "Could not add vector variable '" + name + "' since it already exists.";
            throw std::runtime_error(error_message);
        }

        VectorX vector(rows);
        for (int row = 0; row < vector.rows(); row++)
        {
            vector(row) = Variable(name, row);
        }
        vector_variables.emplace(name, vector);

        return vector;
    }

    MatrixX OptimizationProblem::addVariable(const std::string &name,
                                             size_t rows,
                                             size_t cols)
    {
        if (matrix_variables.find(name) != matrix_variables.end())
        {
            const std::string error_message = "Could not add matrix variable '" + name + "' since it already exists.";
            throw std::runtime_error(error_message);
        }

        MatrixX matrix(rows, cols);
        for (int row = 0; row < matrix.rows(); row++)
        {
            for (int col = 0; col < matrix.cols(); col++)
            {
                matrix(row, col) = Variable(name, row, col);
            }
        }
        matrix_variables.emplace(name, matrix);

        return matrix;
    }

    void OptimizationProblem::addConstraint(const Constraint &constraint)
    {
        if (constraint.getType() == Constraint::Type::Equality)
        {
            this->equality_constraints.push_back(std::get<Constraint::Type::Equality>(constraint.data));
        }
        else if (constraint.getType() == Constraint::Type::Positive)
        {
            this->positive_constraints.push_back(std::get<Constraint::Type::Positive>(constraint.data));
        }
        else if (constraint.getType() == Constraint::Type::Box)
        {
            this->box_constraints.push_back(std::get<Constraint::Type::Box>(constraint.data));
        }
        else if (constraint.getType() == Constraint::Type::SecondOrderCone)
        {
            this->second_order_cone_constraints.push_back(std::get<Constraint::Type::SecondOrderCone>(constraint.data));
        }
    }

    void OptimizationProblem::addConstraint(const std::vector<Constraint> &constraints)
    {
        for (const Constraint &constraint : constraints)
        {
            addConstraint(constraint);
        }
    }

    void OptimizationProblem::addCostTerm(const Scalar &term)
    {
        this->costFunction += term;
    }

    void OptimizationProblem::getVariableValue(const std::string &name, double &var)
    {
        auto found = scalar_variables.find(name);

        if (found != scalar_variables.end())
        {
            var = eval(found->second);
        }
        else
        {
            const std::string error_message = "Could not find scalar variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    void OptimizationProblem::getVariableValue(const std::string &name, Eigen::VectorXd &var)
    {
        auto found = vector_variables.find(name);

        if (found != vector_variables.end())
        {
            var = eval(found->second);
        }
        else
        {
            const std::string error_message = "Could not find vector variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    void OptimizationProblem::getVariableValue(const std::string &name, Eigen::MatrixXd &var)
    {
        auto found = matrix_variables.find(name);

        if (found != matrix_variables.end())
        {
            var = eval(found->second);
        }
        else
        {
            const std::string error_message = "Could not find matrix variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    void OptimizationProblem::getVariable(const std::string &name, Scalar &var)
    {
        auto found = scalar_variables.find(name);

        if (found != scalar_variables.end())
        {
            var = found->second;
        }
        else
        {
            const std::string error_message = "Could not find scalar variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    void OptimizationProblem::getVariable(const std::string &name, VectorX &var)
    {
        auto found = vector_variables.find(name);

        if (found != vector_variables.end())
        {
            var = found->second;
        }
        else
        {
            const std::string error_message = "Could not find vector variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    void OptimizationProblem::getVariable(const std::string &name, MatrixX &var)
    {
        auto found = matrix_variables.find(name);

        if (found != matrix_variables.end())
        {
            var = found->second;
        }
        else
        {
            const std::string error_message = "Could not find matrix variable '" + name + "'. Make sure it has been created first.";
            throw std::runtime_error(error_message);
        }
    }

    double OptimizationProblem::getOptimalValue() const
    {
        return eval(costFunction);
    }

    std::ostream &operator<<(std::ostream &os, const OptimizationProblem &op)
    {
        os << "Minimize\n";
        os << op.costFunction << "\n\n";

        os << "Subject to\n\n";

        os << "Equality Constraints:\n";
        for (const internal::EqualityConstraint &c : op.equality_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Positive Constraints:\n";
        for (const internal::PositiveConstraint &c : op.positive_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Box Constraints:\n";
        for (const internal::BoxConstraint &c : op.box_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Second Order Cone Constraints:\n";
        for (const internal::SecondOrderConeConstraint &c : op.second_order_cone_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";

        return os;
    }

} // namespace cvx