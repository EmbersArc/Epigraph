#pragma once

#include "problem.hpp"

#include <Eigen/Sparse>

namespace cvx
{

    class QPWrapperBase
    {

    public:
        explicit QPWrapperBase(OptimizationProblem &problem);
        virtual bool solve(bool verbose = false) = 0;
        virtual std::string getResultString() const = 0;

        bool isConvex() const;

        size_t getNumVariables() const;
        size_t getNumInequalityConstraints() const;

        friend std::ostream &operator<<(std::ostream &os, const QPWrapperBase &wrapper);

    protected:
        size_t n_variables = 0;

        std::vector<double> solution;

        using MatrixXp = Eigen::Matrix<Parameter, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorXp = Eigen::Matrix<Parameter, Eigen::Dynamic, 1>;

        Eigen::SparseMatrix<Parameter> A_params;
        Eigen::SparseMatrix<Parameter> P_params;
        VectorXp q_params;
        VectorXp l_params;
        VectorXp u_params;

    private:
        void addVariable(Variable &variable);
    };

} // namespace cvx