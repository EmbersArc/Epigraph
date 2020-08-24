/**
 * @file constraint.hpp
 * 
 */

#pragma once

#include "expressions.hpp"

#include <variant>

namespace cvx
{
    namespace internal
    {

        struct EqualityConstraint
        {
            Affine affine;
            friend std::ostream &operator<<(std::ostream &os, const EqualityConstraint &constraint);
        };

        struct PositiveConstraint
        {
            Affine affine;
            friend std::ostream &operator<<(std::ostream &os, const PositiveConstraint &constraint);
        };

        struct BoxConstraint
        {
            Affine lower;
            Affine middle;
            Affine upper;
            friend std::ostream &operator<<(std::ostream &os, const BoxConstraint &constraint);
        };

        struct SecondOrderConeConstraint
        {
            std::vector<Affine> norm;
            Affine affine;
            friend std::ostream &operator<<(std::ostream &os, const SecondOrderConeConstraint &constraint);
        };

    } // namespace internal

    class Constraint
    {
    public:
        enum Type
        {
            Equality,
            Positive,
            Box,
            SecondOrderCone
        };

        Type getType() const;

        friend std::ostream &operator<<(std::ostream &os, const Constraint &constraint);

        friend OptimizationProblem;

        friend Constraint equalTo(const Scalar &lhs, const Scalar &rhs);
        friend Constraint lessThan(const Scalar &lhs, const Scalar &rhs);
        friend Constraint greaterThan(const Scalar &lhs, const Scalar &rhs);
        friend Constraint box(const Scalar &lower, const Scalar &middle, const Scalar &upper);

    private:
        void asEquality(const internal::Affine &affine);
        void asPositive(const internal::Affine &affine);
        void asBox(const internal::Affine &lower, const internal::Affine &middle, const internal::Affine &upper);
        void asSecondOrderCone(const std::vector<internal::Affine> &norm, const internal::Affine &affine);

        using constraint_variant_t = std::variant<internal::EqualityConstraint,
                                                  internal::PositiveConstraint,
                                                  internal::BoxConstraint,
                                                  internal::SecondOrderConeConstraint>;
        constraint_variant_t data;
    };

    /**
     * @brief Create an equality constraint: lhs == rhs
     * 
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    Constraint equalTo(const Scalar &lhs, const Scalar &rhs);

    /**
     * @brief Create a greater than or equal constraint: lhs >= rhs
     * 
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    Constraint greaterThan(const Scalar &lhs, const Scalar &rhs);

    /**
     * @brief Create a less than or equal constraint: lhs <= rhs
     * 
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    Constraint lessThan(const Scalar &lhs, const Scalar &rhs);

    /**
     * @brief Create a box constraint: lower <= middle <= upper 
     * 
     * @param lower  The lower bound
     * @param middle The middle term
     * @param upper  The upper bound
     * @return Constraint The constraint to pass to addConstraint()
     */
    Constraint box(const Scalar &lower, const Scalar &middle, const Scalar &upper);

    template <typename Derived>
    std::vector<Constraint> equalTo(const Eigen::MatrixBase<Derived> &lhs, const Scalar &rhs)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<Derived>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        for (int row = 0; row < lhs.rows(); row++)
        {
            for (int col = 0; col < lhs.cols(); col++)
            {
                constraints.push_back(equalTo(lhs(row, col), rhs));
            }
        }

        return constraints;
    }

    /**
     * @brief Create an equality constraint: lhs == rhs
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> equalTo(const Scalar &lhs, const Eigen::MatrixBase<Derived> &rhs)
    {
        return equalTo(rhs, lhs);
    }

    /**
     * @brief Create an equality constraint: lhs == rhs
     * 
     * @tparam DerivedLhs Has to be cvx::Scalar
     * @tparam DerivedRhs Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename DerivedLhs, typename DerivedRhs>
    std::vector<Constraint> equalTo(const Eigen::MatrixBase<DerivedLhs> &lhs, const Eigen::MatrixBase<DerivedRhs> &rhs)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedLhs>::Scalar, Scalar>);
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedRhs>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        // If they are both matrices then the dimensions have to match
        if (lhs.rows() != rhs.rows() or lhs.cols() != rhs.cols())
        {
            throw std::runtime_error("Invalid dimensions in constraint.");
        }

        for (int row = 0; row < rhs.rows(); row++)
        {
            for (int col = 0; col < rhs.cols(); col++)
            {
                constraints.push_back(equalTo(lhs(row, col), rhs(row, col)));
            }
        }

        return constraints;
    }

    /**
     * @brief Create a less than or equal constraint: lhs <= rhs
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> lessThan(const Eigen::MatrixBase<Derived> &lhs, const Scalar &rhs)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<Derived>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        for (int row = 0; row < lhs.rows(); row++)
        {
            for (int col = 0; col < lhs.cols(); col++)
            {
                constraints.push_back(lessThan(lhs(row, col), rhs));
            }
        }

        return constraints;
    }

    /**
     * @brief Create a less than or equal constraint: lhs <= rhs
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> lessThan(const Scalar &lhs, const Eigen::MatrixBase<Derived> &rhs)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<Derived>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        for (int row = 0; row < rhs.rows(); row++)
        {
            for (int col = 0; col < rhs.cols(); col++)
            {
                constraints.push_back(lessThan(lhs, rhs(row, col)));
            }
        }

        return constraints;
    }

    /**
     * @brief Create a less than or equal constraint: lhs <= rhs
     * 
     * @tparam DerivedLhs Has to be cvx::Scalar
     * @tparam DerivedRhs Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename DerivedLhs, typename DerivedRhs>
    std::vector<Constraint> lessThan(const Eigen::MatrixBase<DerivedLhs> &lhs, const Eigen::MatrixBase<DerivedRhs> &rhs)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedLhs>::Scalar, Scalar>);
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedRhs>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        if ((lhs.rows() != rhs.rows()) or
            (lhs.cols() != rhs.cols()))
        {
            throw std::runtime_error("Invalid dimensions in constraint.");
        }

        for (int row = 0; row < rhs.rows(); row++)
        {
            for (int col = 0; col < rhs.cols(); col++)
            {
                constraints.push_back(lessThan(lhs(row, col), rhs(row, col)));
            }
        }

        return constraints;
    }

    /**
     * @brief Create a greater than or equal constraint: lhs >= rhs
     * 
     * @tparam DerivedLhs Has to be cvx::Scalar
     * @tparam DerivedRhs Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename DerivedLhs, typename DerivedRhs>
    std::vector<Constraint> greaterThan(const Eigen::MatrixBase<DerivedLhs> &lhs, const Eigen::MatrixBase<DerivedRhs> &rhs)
    {
        return lessThan(rhs, lhs);
    }

    /**
     * @brief Create a greater than or equal constraint: lhs >= rhs
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> greaterThan(const Scalar &lhs, const Eigen::MatrixBase<Derived> &rhs)
    {
        return lessThan(rhs, lhs);
    }

    /**
     * @brief Create a greater than or equal constraint: lhs >= rhs
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lhs The left hand side
     * @param rhs The right hand side
     * @return Constraint The constraint to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> greaterThan(const Eigen::MatrixBase<Derived> &lhs, const Scalar &rhs)
    {
        return lessThan(rhs, lhs);
    }

    /**
     * @brief Create a box constraint: lower <= middle <= upper 
     * 
     * @tparam DerivedLower Has to be cvx::Scalar
     * @tparam DerivedMiddle Has to be cvx::Scalar
     * @tparam DerivedUpper Has to be cvx::Scalar
     * @param lower  The lower bound
     * @param middle The middle term
     * @param upper  The upper bound
     * @return std::vector<Constraint> 
     */
    template <typename DerivedLower, typename DerivedMiddle, typename DerivedUpper>
    std::vector<Constraint> box(const Eigen::MatrixBase<DerivedLower> &lower,
                                const Eigen::MatrixBase<DerivedMiddle> &middle,
                                const Eigen::MatrixBase<DerivedUpper> &upper)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedLower>::Scalar, Scalar>);
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedMiddle>::Scalar, Scalar>);
        static_assert(std::is_same_v<typename Eigen::MatrixBase<DerivedUpper>::Scalar, Scalar>);

        if (lower.rows() != middle.rows() or middle.rows() != upper.rows() or
            lower.cols() != middle.cols() or middle.cols() != upper.cols())
        {
            throw std::runtime_error("Invalid dimensions in constraint.");
        }

        std::vector<Constraint> constraints;

        for (int row = 0; row < middle.rows(); row++)
        {
            for (int col = 0; col < middle.cols(); col++)
            {
                constraints.push_back(box(lower(row, col),
                                          middle(row, col),
                                          upper(row, col)));
            }
        }

        return constraints;
    }

    /**
     * @brief Create a box constraint: lower <= middle <= upper 
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param lower  The lower bound
     * @param middle The middle term
     * @param upper  The upper bound
     * @return std::vector<Constraint> The constraints to pass to addConstraint()
     */
    template <typename Derived>
    std::vector<Constraint> box(const Scalar &lower,
                                const Eigen::MatrixBase<Derived> &middle,
                                const Scalar &upper)
    {
        static_assert(std::is_same_v<typename Eigen::MatrixBase<Derived>::Scalar, Scalar>);

        std::vector<Constraint> constraints;

        for (int row = 0; row < middle.rows(); row++)
        {
            for (int col = 0; col < middle.cols(); col++)
            {
                constraints.push_back(box(lower,
                                          middle(row, col),
                                          upper));
            }
        }

        return constraints;
    }

} // namespace cvx