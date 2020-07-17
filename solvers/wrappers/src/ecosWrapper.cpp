#include "ecosWrapper.hpp"

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
            current_G->data(),
            G_col_ind.data(),
            G_row_ind.data(),
            current_A->data(),
            A_col_ind.data(),
            A_row_ind.data(),
            current_c->data(),
            current_h->data(),
            current_b->data());

        if (work == nullptr)
        {
            throw std::runtime_error("ECOS failed to set up problem.");
        }
    }

    bool ECOSSolver::solve(bool verbose)
    {
        work->stgs->verbose = verbose;

        update();

        ECOS_updateData(work,
                        current_G->data(),
                        current_A->data(),
                        current_c->data(),
                        current_h->data(),
                        current_b->data());

        exitflag = ECOS_solve(work);

        solution = std::vector<double>(work->x, work->x + getNumVariables());

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

    ECOSSolver::~ECOSSolver()
    {
        ECOS_cleanup(work, 0);
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
        // Awkward switching between memory locations
        current_G = alternate_memory ? &G1 : &G2;
        current_A = alternate_memory ? &A1 : &A2;
        current_c = alternate_memory ? &c1 : &c2;
        current_h = alternate_memory ? &h1 : &h2;
        current_b = alternate_memory ? &b1 : &b2;

        // The signs for A and G must be flipped because they are negative in the ECOS interface
        *current_G = -eval(Eigen::Map<VectorXp>(G_params.valuePtr(), G_params.nonZeros()));
        *current_A = -eval(Eigen::Map<VectorXp>(A_params.valuePtr(), A_params.nonZeros()));
        *current_c = eval(c_params);
        *current_h = eval(h_params);
        *current_b = eval(b_params);

        alternate_memory = !alternate_memory;
    }

} // namespace cvx::ecos
