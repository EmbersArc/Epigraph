
#pragma once

#include "problem.hpp"

#include <Eigen/Sparse>

namespace cvx::internal
{

    class WrapperBase
    {

    public:
        virtual bool solve(bool verbose = false) = 0;
        virtual std::string getResultString() const = 0;
        size_t getNumVariables() const;

        virtual ~WrapperBase();

    protected:
        using MatrixXp = Eigen::Matrix<Parameter, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorXp = Eigen::Matrix<Parameter, Eigen::Dynamic, 1>;
        std::vector<Variable> variables;
        std::vector<double> solution;

        virtual void addVariable(Variable &variable) = 0;

    private:
    };

} // namespace cvx