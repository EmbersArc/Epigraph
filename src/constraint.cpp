#include "constraint.hpp"

namespace cvx
{

    std::ostream &operator<<(std::ostream &os, const EqualityConstraint &constraint)
    {
        os << constraint.affine << " == 0";
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const PositiveConstraint &constraint)
    {
        os << "0 <= " << constraint.affine;
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const BoxConstraint &constraint)
    {
        os << constraint.lower << " <= " << constraint.middle << " <= " << constraint.upper;
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const SecondOrderConeConstraint &constraint)
    {
        os << "(";
        for (size_t i = 0; i < constraint.norm.size(); i++)
        {
            os << "(" << constraint.norm[i] << ")^2";
            if (i != constraint.norm.size() - 1)
            {
                os << " + ";
            }
        }
        os << ")^(1/2)";

        os << " <= " << constraint.affine;

        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Constraint &constraint)
    {
        switch (constraint.data.index())
        {
        case Constraint::Type::Equality:
            os << std::get<Constraint::Type::Equality>(constraint.data);
            break;
        case Constraint::Type::Positive:
            os << std::get<Constraint::Type::Positive>(constraint.data);
            break;
        case Constraint::Type::Box:
            os << std::get<Constraint::Type::Box>(constraint.data);
            break;
        case Constraint::Type::SecondOrderCone:
            os << std::get<Constraint::Type::SecondOrderCone>(constraint.data);
            break;
        }
        return os;
    }

    Constraint::Type Constraint::getType() const
    {
        return Type(data.index());
    }

    // affine == 0
    void Constraint::asEquality(const Affine &affine)
    {
        EqualityConstraint constraint;
        constraint.affine = affine;
        data = constraint;
    }

    // affine >= 0
    void Constraint::asPositive(const Affine &affine)
    {
        PositiveConstraint constraint;
        constraint.affine = affine;
        data = constraint;
    }

    // affine <= affine <= affine
    void Constraint::asBox(const Affine &lower, const Affine &middle, const Affine &upper)
    {
        BoxConstraint constraint;
        constraint.lower = lower;
        constraint.middle = middle;
        constraint.upper = upper;
        data = constraint;
    }

    // norm <= affine
    void Constraint::asSecondOrderCone(const std::vector<Affine> &norm, const Affine &affine)
    {
        SecondOrderConeConstraint constraint;
        constraint.norm = norm;
        constraint.affine = affine;
        data = constraint;
    }

    Constraint equalTo(const Scalar &lhs, const Scalar &rhs)
    {
        if (lhs.getOrder() > 1 or rhs.getOrder() > 1)
        {
            throw std::runtime_error("The terms in an equality have to be constant or linear.");
        }
        Constraint constraint;
        constraint.asEquality(lhs.affine - rhs.affine);
        return constraint;
    }

    Constraint lessThan(const Scalar &lhs, const Scalar &rhs)
    {
        if (rhs.getOrder() > 1)
        {
            throw std::runtime_error("The larger term in an inequality has to be constant or linear.");
        }

        Constraint constraint;
        if (lhs.isNorm())
        {
            std::vector<Affine> norm2_terms;
            for (const Product &product : lhs.products)
            {
                norm2_terms.push_back(product.firstTerm());
            }
            constraint.asSecondOrderCone(norm2_terms, rhs.affine - lhs.affine);
        }
        else if (lhs.getOrder() < 2)
        {
            if (lhs.getOrder() > 0 or rhs.getOrder() > 0)
                constraint.asPositive(rhs.affine - lhs.affine);
        }
        else
        {
            throw std::runtime_error("The smaller term in an inequality has to be constant, linear or a 2-norm.");
        }

        return constraint;
    }

    Constraint greaterThan(const Scalar &lhs, const Scalar &rhs)
    {
        return lessThan(rhs, lhs);
    }

    Constraint box(const Scalar &lower, const Scalar &middle, const Scalar &upper)
    {
        if (lower.getOrder() > 1 or middle.getOrder() > 1 or upper.getOrder() > 1)
        {
            throw std::runtime_error("The terms in box constraints have to be constant or linear.");
        }

        Constraint constraint;
        constraint.asBox(lower.affine, middle.affine, upper.affine);
        return constraint;
    }

} // namespace cvx
