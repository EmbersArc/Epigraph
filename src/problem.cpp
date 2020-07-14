#include "problem.hpp"

namespace cvx
{
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
        for (const EqualityConstraint &c : op.equality_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Positive Constraints:\n";
        for (const PositiveConstraint &c : op.positive_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Box Constraints:\n";
        for (const BoxConstraint &c : op.box_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";
        os << "Second Order Cone Constraints:\n";
        for (const SecondOrderConeConstraint &c : op.second_order_cone_constraints)
        {
            os << c << "\n\n";
        }
        os << "\n";

        return os;
    }

} // namespace cvx