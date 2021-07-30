#pragma once

#include "problem.hpp"

namespace cvx::internal
{

    class WrapperBase
    {

    public:
        WrapperBase() = default;
        virtual ~WrapperBase();
        virtual bool solve(bool verbose = false) = 0;
        virtual std::string getResultString() const = 0;
        size_t getNumVariables() const;
        virtual bool isFeasible(double tolerance) const = 0;

    protected:
        using MatrixXp = Eigen::Matrix<Parameter, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorXp = Eigen::Matrix<Parameter, Eigen::Dynamic, 1>;
        std::vector<Variable> variables;
        std::shared_ptr<std::vector<double>> solution = std::make_shared<std::vector<double>>();

        virtual void addVariable(Variable &variable) = 0;

    private:
        WrapperBase(const WrapperBase &);
        static size_t solver_count;
    };

} // namespace cvx::internal