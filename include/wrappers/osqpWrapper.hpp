#pragma once

#include "wrappers/qpWrapperBase.hpp"

#include "osqp.h"

namespace cvx::osqp
{

    class OSQPSolver final : public internal::QPWrapperBase
    {

    public:
        explicit OSQPSolver(OptimizationProblem &problem);
        bool solve(bool verbose = false) override;
        std::string getResultString() const override;
        const OSQPSettings &getSettings() const;
        const OSQPInfo &getInfo() const;
        c_int getExitCode() const;

        // settings
        void setAlpha(c_float alpha);
        void setDelta(c_float delta);
        void setEpsAbs(c_float eps);
        void setEpsPinf(c_float eps);
        void setEpsDinf(c_float eps);
        void setEpsRel(c_float eps);
        void setMaxIter(c_int iter);
        void setPolish(bool polish);
        void setPolishRefine(c_int iter);
        void setRho(c_float rho);
        void setScaledTermination(bool scaled_termination);
        void setTimeLimit(c_float t);
        void setCheckTermination(c_int interval);
        void setWarmStart(bool warm_start);

        ~OSQPSolver();

    private:
        Eigen::SparseMatrix<c_float> P;
        Eigen::SparseMatrix<c_float> A;
        Eigen::Matrix<c_float, Eigen::Dynamic, 1> q;
        Eigen::Matrix<c_float, Eigen::Dynamic, 1> l;
        Eigen::Matrix<c_float, Eigen::Dynamic, 1> u;

        c_int exitflag = OSQP_UNSOLVED;

        Eigen::Matrix<c_int, Eigen::Dynamic, 1> P_row_ind;
        Eigen::Matrix<c_int, Eigen::Dynamic, 1> P_col_ind;
        Eigen::Matrix<c_int, Eigen::Dynamic, 1> A_row_ind;
        Eigen::Matrix<c_int, Eigen::Dynamic, 1> A_col_ind;

        OSQPWorkspace *workspace;
        OSQPData data;
        OSQPSettings settings;

        void update();
        void cleanUp();
    };

} // namespace cvx::osqp