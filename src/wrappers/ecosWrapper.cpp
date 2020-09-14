#include "wrappers/ecosWrapper.hpp"

namespace cvx::ecos
{

    ECOSSolver::ECOSSolver(OptimizationProblem &problem) : SOCPWrapperBase(problem)
    {
        update();

        cone_constraint_dimensions = soc_dims.cast<idxint>();
        G_row_ind = Eigen::Map<Eigen::VectorXi>(G_params.innerIndexPtr(), G_params.nonZeros()).cast<idxint>();
        A_row_ind = Eigen::Map<Eigen::VectorXi>(A_params.innerIndexPtr(), A_params.nonZeros()).cast<idxint>();
        G_col_ind = Eigen::Map<Eigen::VectorXi>(G_params.outerIndexPtr(), G_params.outerSize() + 1).cast<idxint>();
        A_col_ind = Eigen::Map<Eigen::VectorXi>(A_params.outerIndexPtr(), A_params.outerSize() + 1).cast<idxint>();

        work = ECOS_setup(
            getNumVariables(),
            getNumInequalityConstraints(),
            getNumEqualityConstraints(),
            getNumPositiveConstraints(),
            getNumCones(),
            cone_constraint_dimensions.data(),
            0,
            G.data(),
            G_col_ind.data(),
            G_row_ind.data(),
            A.data(),
            A_col_ind.data(),
            A_row_ind.data(),
            c.data(),
            h.data(),
            b.data());

        if (work == nullptr)
        {
            cleanUp();
            throw std::runtime_error("ECOS failed to set up problem.");
        }
    }

    bool ECOSSolver::solve(bool verbose)
    {
        work->stgs->verbose = verbose;

        update();

        ECOS_updateData(work,
                        G.data(),
                        A.data(),
                        c.data(),
                        h.data(),
                        b.data());

        exitflag = ECOS_solve(work);

        *solution = std::vector<double>(work->x, work->x + getNumVariables());

        if (exitflag == ECOS_SIGINT)
        {
            std::exit(130);
        }

        return exitflag != ECOS_FATAL;
    }

    std::string ECOSSolver::getResultString() const
    {
        switch (exitflag)
        {
        case ECOS_UNSOLVED:
            return "Problem not solved yet.";
        case ECOS_OPTIMAL:
            return "Optimal solution found.";
        case ECOS_PINF:
            return "Certificate of primal infeasibility found.";
        case ECOS_DINF:
            return "Certificate of dual infeasibility found.";

        case ECOS_OPTIMAL + ECOS_INACC_OFFSET:
            return "Optimal solution found subject to reduced tolerances.";
        case ECOS_PINF + ECOS_INACC_OFFSET:
            return "Certificate of primal infeasibility found subject to reduced tolerances.";
        case ECOS_DINF + ECOS_INACC_OFFSET:
            return "Certificate of dual infeasibility found subject to reduced tolerances.";

        case ECOS_MAXIT:
            return "Maximum number of iterations reached.";

        case ECOS_NUMERICS:
            return "Numerical problems (unreliable search direction).";
        case ECOS_OUTCONE:
            return "Numerical problems (slacks or multipliers outside cone)";

        case ECOS_SIGINT:
            return "Interrupted by signal or CTRL-C.";
        default: // ECOS_FATAL:
            return "Unknown problem in solver.";
        }
    }

    settings &ECOSSolver::getSettings()
    {
        return *work->stgs;
    }

    const stats &ECOSSolver::getInfo() const
    {
        return *work->info;
    }

    idxint ECOSSolver::getExitCode() const
    {
        return exitflag;
    }

    void ECOSSolver::update()
    {
        // The signs for A and G must be flipped because they are negative in the ECOS interface
        G = -eval(Eigen::Map<VectorXp>(G_params.valuePtr(), G_params.nonZeros()));
        A = -eval(Eigen::Map<VectorXp>(A_params.valuePtr(), A_params.nonZeros()));
        c = eval(c_params);
        h = eval(h_params);
        b = eval(b_params);
    }

    void ECOSSolver::cleanUp()
    {
        ECOS_cleanup(work, 0);
    }

    ECOSSolver::~ECOSSolver()
    {
        cleanUp();
    }

} // namespace cvx::ecos
