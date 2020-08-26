#pragma once

#include "socpWrapperBase.hpp"

#define DCTRLC = 1
#define DLONG
#define LDL_LONG
#include "ecos/include/ecos.h"
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

        bool alternate_memory = false;

        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> G1;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> A1;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> c1;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> h1;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> b1;

        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> G2;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> A2;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> c2;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> h2;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> b2;

        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> *current_G;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> *current_A;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> *current_c;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> *current_h;
        Eigen::Matrix<pfloat, Eigen::Dynamic, 1> *current_b;

        Eigen::Matrix<idxint, Eigen::Dynamic, 1> cone_constraint_dimensions;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> G_row_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> G_col_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> A_row_ind;
        Eigen::Matrix<idxint, Eigen::Dynamic, 1> A_col_ind;
    };

} // namespace cvx::ecos