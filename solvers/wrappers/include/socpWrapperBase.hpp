#pragma once

#include "wrapperBase.hpp"

namespace cvx::internal
{

    class SOCPWrapperBase : public WrapperBase
    {

    public:
        explicit SOCPWrapperBase(OptimizationProblem &problem);

        size_t getNumEqualityConstraints() const;
        size_t getNumInequalityConstraints() const;
        size_t getNumPositiveConstraints() const;
        size_t getNumCones() const;

        friend std::ostream &operator<<(std::ostream &os, const SOCPWrapperBase &wrapper);

    protected:
        Eigen::SparseMatrix<Parameter> A_params;
        Eigen::SparseMatrix<Parameter> G_params;
        VectorXp c_params;
        VectorXp h_params;
        VectorXp b_params;
        Eigen::VectorXi soc_dims;

    private:
        void addVariable(Variable &variable) final override;
    };

} // namespace cvx