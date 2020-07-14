#pragma once

#include "socpWrapperBase.hpp"

#include "eicos/include/eicos.hpp"

namespace cvx::eicos
{

    class EiCOSSolver final : public SOCPWrapperBase
    {
    public:
        explicit EiCOSSolver(OptimizationProblem &problem);
        bool solve(bool verbose = false) override;
        std::string getResultString() const override;
        EiCOS::Settings &getSettings();
        const EiCOS::Information &getInfo() const;
        EiCOS::exitcode getExitCode() const;

    private:
        EiCOS::exitcode exitflag;

        std::unique_ptr<EiCOS::Solver> solver;

        Eigen::SparseMatrix<double> A;
        Eigen::SparseMatrix<double> G;
        Eigen::VectorXd c;
        Eigen::VectorXd h;
        Eigen::VectorXd b;

        void update();
    };

} // namespace cvx