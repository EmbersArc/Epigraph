#pragma once

#include "wrappers/socpWrapperBase.hpp"

#include "ecos.h"
#define ECOS_UNSOLVED (-10)

namespace cvx::ecos
{
    class ECOSSolver final : public internal::SOCPWrapperBase
    {

    public:
        explicit ECOSSolver(OptimizationProblem &problem);
        bool solve(bool verbose = false) override;
        std::string getResultString() const override;
        settings &getSettings();
        const stats &getInfo() const;
        idxint getExitCode() const;
        ~ECOSSolver();

    private:
        void update();
        void cleanUp();

        idxint exitflag = ECOS_UNSOLVED;

        pwork *work;

        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> G;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> A;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> c;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> h;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> b;

        Eigen::Matrix<idxint, Eigen::Dynamic, 1> cone_constraint_dimensions;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> G_row_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> G_col_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> A_row_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> A_col_ind;
    };

} // namespace cvx::ecos