#pragma once

#include "problem.hpp"

#include <Eigen/Sparse>

namespace cvx
{

    class SOCPWrapperBase
    {

    public:
        explicit SOCPWrapperBase(OptimizationProblem &problem);
        virtual bool solve(bool verbose = false) = 0;
        virtual std::string getResultString() const = 0;

        size_t getNumEqualityConstraints() const;
        size_t getNumInequalityConstraints() const;
        size_t getNumPositiveConstraints() const;
        size_t getNumCones() const;

        friend std::ostream &operator<<(std::ostream &os, const SOCPWrapperBase &wrapper);

    protected:
        size_t n_variables = 0;

        std::vector<double> solution;

        using MatrixXp = Eigen::Matrix<Parameter, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorXp = Eigen::Matrix<Parameter, Eigen::Dynamic, 1>;
        Eigen::SparseMatrix<Parameter> A_params;
        Eigen::SparseMatrix<Parameter> G_params;
        VectorXp c_params;
        VectorXp h_params;
        VectorXp b_params;
        Eigen::VectorXi soc_dims;
        size_t getNumVariables() const;

    private:
        void addVariable(Variable &variable);
    };

} // namespace cvx