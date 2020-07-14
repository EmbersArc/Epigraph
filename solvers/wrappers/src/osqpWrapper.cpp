#include "osqpWrapper.hpp"

namespace cvx::osqp
{

    OSQPSolver::OSQPSolver(OptimizationProblem &problem) : QPWrapperBase(problem)
    {
        update();

        P_row_ind = Eigen::Map<Eigen::VectorXi>(P.innerIndexPtr(), P.nonZeros()).cast<c_int>();
        A_row_ind = Eigen::Map<Eigen::VectorXi>(A.innerIndexPtr(), A.nonZeros()).cast<c_int>();
        P_col_ind = Eigen::Map<Eigen::VectorXi>(P.outerIndexPtr(), P.outerSize() + 1).cast<c_int>();
        A_col_ind = Eigen::Map<Eigen::VectorXi>(A.outerIndexPtr(), A.outerSize() + 1).cast<c_int>();

        // Populate data
        data.n = getNumVariables();
        data.m = getNumInequalityConstraints();
        data.P = csc_matrix(getNumVariables(), getNumVariables(), P.nonZeros(),
                            P.valuePtr(), P_row_ind.data(), P_col_ind.data());
        data.q = q.data();
        data.A = csc_matrix(getNumInequalityConstraints(), getNumVariables(), A.nonZeros(),
                            A.valuePtr(), A_row_ind.data(), A_col_ind.data());
        data.l = l.data();
        data.u = u.data();

        osqp_set_default_settings(&settings);

        osqp_setup(&workspace, &data, &settings);
    }

    void OSQPSolver::update()
    {
        P = eval(P_params);
        q = eval(q_params);
        A = eval(A_params);
        l = eval(l_params);
        u = eval(u_params);
    }

    bool OSQPSolver::solve(bool verbose)
    {
        osqp_update_verbose(workspace, verbose);

        update();

        osqp_update_P(workspace, P.valuePtr(), OSQP_NULL, 0);
        osqp_update_A(workspace, A.valuePtr(), OSQP_NULL, 0);
        osqp_update_lin_cost(workspace, q.data());
        osqp_update_bounds(workspace, l.data(), u.data());

        exitflag = osqp_solve(workspace);

        solution = std::vector<double>(workspace->solution->x,
                                       workspace->solution->x + getNumVariables());

        return exitflag == 0;
    }

    std::string OSQPSolver::getResultString() const
    {
        return workspace->info->status;
    }

    const OSQPSettings &OSQPSolver::getSettings() const
    {
        return settings;
    }

    const OSQPInfo &OSQPSolver::getInfo() const
    {
        return *workspace->info;
    }

    osqp::c_int OSQPSolver::getExitCode() const
    {
        return exitflag;
    }

    void OSQPSolver::setAlpha(c_float alpha)
    {
        osqp_update_alpha(workspace, alpha);
    }

    void OSQPSolver::setDelta(c_float delta)
    {
        osqp_update_delta(workspace, delta);
    }

    void OSQPSolver::setEpsAbs(c_float eps)
    {
        osqp_update_eps_abs(workspace, eps);
    }

    void OSQPSolver::setEpsPinf(c_float eps)
    {
        osqp_update_eps_prim_inf(workspace, eps);
    }

    void OSQPSolver::setEpsDinf(c_float eps)
    {
        osqp_update_eps_dual_inf(workspace, eps);
    }

    void OSQPSolver::setEpsRel(c_float eps)
    {
        osqp_update_eps_rel(workspace, eps);
    }

    void OSQPSolver::setMaxIter(c_int iter)
    {
        osqp_update_max_iter(workspace, iter);
    }

    void OSQPSolver::setPolish(bool polish)
    {
        osqp_update_polish(workspace, c_int(polish));
    }

    void OSQPSolver::setPolishRefine(c_int iter)
    {
        osqp_update_polish_refine_iter(workspace, iter);
    }

    void OSQPSolver::setRho(c_float rho)
    {
        osqp_update_rho(workspace, rho);
    }

    void OSQPSolver::setScaledTermination(bool scaled_termination)
    {
        osqp_update_scaled_termination(workspace, c_int(scaled_termination));
    }

    void OSQPSolver::setTimeLimit(c_float t)
    {
        osqp_update_time_limit(workspace, t);
    }

    void OSQPSolver::setCheckTermination(c_int interval)
    {
        osqp_update_check_termination(workspace, interval);
    }

    void OSQPSolver::setWarmStart(bool warm_start)
    {
        osqp_update_warm_start(workspace, c_int(warm_start));
    }

    OSQPSolver::~OSQPSolver()
    {
        osqp_cleanup(workspace);
    }

} // namespace cvx::osqp
