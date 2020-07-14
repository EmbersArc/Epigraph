#pragma once

#include "socpWrapperBase.hpp"

namespace cvx::ecos
{
#define DCTRLC = 1
#define DLONG
#define LDL_LONG
#include "ecos/include/ecos.h"
#define ECOS_UNSOLVED (-10)
    static_assert(std::is_same_v<idxint, SuiteSparse_long> and std::is_same_v<idxint, long>,
                  "Definitions of idxint are not consistent."
                  "Make sure ECOS is compiled with USE_LONG = 1.");

    class ECOSSolver final : public SOCPWrapperBase
    {

    public:
        explicit ECOSSolver(OptimizationProblem &problem);
        bool solve(bool verbose = false);
        std::string getResultString() const;
        settings &getSettings();
        const stats &getInfo() const;
        idxint getExitCode() const;
        ~ECOSSolver();

    private:
        void update();

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