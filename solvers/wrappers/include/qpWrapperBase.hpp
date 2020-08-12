#pragma once

#include "wrapperBase.hpp"

namespace cvx::internal
{

    class QPWrapperBase : public WrapperBase
    {

    public:
        explicit QPWrapperBase(OptimizationProblem &problem);

        bool isConvex() const;

        size_t getNumInequalityConstraints() const;

        friend std::ostream &operator<<(std::ostream &os, const QPWrapperBase &wrapper);

    protected:
        Eigen::SparseMatrix<Parameter> A_params;
        Eigen::SparseMatrix<Parameter> P_params;
        VectorXp q_params;
        VectorXp l_params;
        VectorXp u_params;

    private:
        void addVariable(Variable &variable) final override;
    };

} // namespace cvx