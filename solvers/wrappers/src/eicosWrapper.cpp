#include "eicosWrapper.hpp"

namespace cvx::eicos
{

    EiCOSSolver::EiCOSSolver(OptimizationProblem &problem) : SOCPWrapperBase(problem)
    {
        update();

        solver = std::make_unique<EiCOS::Solver>(G, A, c, h, b, soc_dims);
    }

    void EiCOSSolver::update()
    {
        // The signs for A and G must be flipped because they are negative in the EiCOS interface
        G = -eval(G_params);
        A = -eval(A_params);
        c = eval(c_params);
        h = eval(h_params);
        b = eval(b_params);
    }

    bool EiCOSSolver::solve(bool verbose)
    {
        update();

        solver->updateData(G, A, c, h, b);

        exitflag = solver->solve(verbose);

        // copy solution
        solution = std::vector<double>(solver->solution().data(),
                                       solver->solution().data() + solver->solution().size());

        return exitflag != EiCOS::exitcode::fatal;
    }

    std::string EiCOSSolver::getResultString() const
    {
        switch (exitflag)
        {
        case EiCOS::exitcode::optimal:
            return "Optimal solution found.";
        case EiCOS::exitcode::primal_infeasible:
            return "Certificate of primal infeasibility found.";
        case EiCOS::exitcode::dual_infeasible:
            return "Certificate of dual infeasibility found.";

        case EiCOS::exitcode::close_to_optimal:
            return "Optimal solution found subject to reduced tolerances.";
        case EiCOS::exitcode::close_to_primal_infeasible:
            return "Certificate of primal infeasibility found subject to reduced tolerances.";
        case EiCOS::exitcode::close_to_dual_infeasible:
            return "Certificate of dual infeasibility found subject to reduced tolerances.";

        case EiCOS::exitcode::maxit:
            return "Maximum number of iterations reached.";

        case EiCOS::exitcode::numerics:
            return "Numerical problems (unreliable search direction).";
        case EiCOS::exitcode::outcone:
            return "Numerical problems (slacks or multipliers outside cone)";

        default: // EiCOS::exitcode::fatal
            return "Unknown problem in solver.";
        }
    }

    EiCOS::Settings &EiCOSSolver::getSettings()
    {
        return solver->getSettings();
    }

    const EiCOS::Information &EiCOSSolver::getInfo() const
    {
        return solver->getInfo();
    }

    EiCOS::exitcode EiCOSSolver::getExitCode() const
    {
        return exitflag;
    }

} // namespace cvx::eicos
