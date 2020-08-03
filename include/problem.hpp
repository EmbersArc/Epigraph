#pragma once

#include "constraint.hpp"

#include <map>

namespace cvx
{

    class SOCPWrapperBase;
    class QPWrapperBase;

    class OptimizationProblem
    {
    public:
        /**
         * @brief Creates and returns a single variable.
         * 
         * @param name The name of the variable
         * @return Scalar The variable
         */
        Scalar addVariable(const std::string &name);

        /**
         * @brief Creates and returns a vector of variables.
         * 
         * @param name The name of the variable
         * @param rows The number of elements in the vector
         * @return VectorX The vector of variables
         */
        VectorX addVariable(const std::string &name,
                            size_t rows);

        /**
         * @brief Creates and returns a matrix of variables.
         * 
         * @param name The name of the variable
         * @param rows The number of rows of the matrix
         * @param cols The number of columns in the matrix
         * @return MatrixX The matrix of variables
         */
        MatrixX addVariable(const std::string &name,
                            size_t rows,
                            size_t cols);

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
         * @brief Get a variable that exists in the problem.
         * 
         * @param name The name of the variable
         * @param var The retured variable
         */
        void getVariable(const std::string &name, Scalar &var);
        void getVariable(const std::string &name, VectorX &var);
        void getVariable(const std::string &name, MatrixX &var);

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
        std::map<std::string, Variable> variables;
    };

} // namespace cvx