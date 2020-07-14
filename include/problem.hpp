#pragma once

#include "constraint.hpp"

namespace cvx
{

    class SOCPWrapperBase;
    class QPWrapperBase;

    class OptimizationProblem
    {
    public:
        /**
         * @brief Add a single constraint to the problem.
         * 
         * @param constraint The constraint created by equalTo(), lessThan(), greaterThan() or box()
         */
        void addConstraint(const Constraint &constraint);

        /**
         * @brief Add multiple constraints to the problem.
         * 
         * @param constraint The constraint created by equalTo(), lessThan(), greaterThan() or box()
         */
        void addConstraint(const std::vector<Constraint> &constraints);

        /**
         * @brief Add a cost term to the problem's cost function. This has to be a scalar.
         * 
         * @param term A scalar cost term
         */
        void addCostTerm(const Scalar &term);

        /**
         * @brief Returns the evaluated cost function. Only call this after the problem has been solved.
         * 
         * @note This value may differ from the optimal value reported by the solver.
         * 
         * @return double 
         */
        double getOptimalValue() const;

        /**
         * @brief Returns the number of used variables in the problem.
         * 
         * @return size_t The number of variables
         */
        size_t getNumVariables() const;

        friend std::ostream &operator<<(std::ostream &os, const OptimizationProblem &socp);
        friend SOCPWrapperBase;
        friend QPWrapperBase;

    private:
        Scalar costFunction;
        std::vector<EqualityConstraint> equality_constraints;
        std::vector<PositiveConstraint> positive_constraints;
        std::vector<BoxConstraint> box_constraints;
        std::vector<SecondOrderConeConstraint> second_order_cone_constraints;
    };

} // namespace cvx