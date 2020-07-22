import cvxpy as cp
import numpy as np

# Generate a random non-trivial quadratic program.
m = 15
n = 10
p = 5
np.random.seed(1)
P = np.random.randn(n, n)
P = P.T @ P
q = np.random.randn(n)
G = np.random.randn(m, n)
h = G @ np.random.randn(n)
A = np.random.randn(p, n)
b = np.random.randn(p)

# Define and solve the CVXPY problem.
x = cp.Variable(n)
prob = cp.Problem(cp.Minimize(cp.quad_form(x, 0.5 * P) + q.T @ x),
                  [G @ x <= h,
                   A @ x == b])
prob.solve(solver=cp.OSQP, verbose=True)

# Print result.
print("\nThe optimal value is", prob.value)
print("A solution x is")
print(repr(x.value))
print("A dual solution corresponding to the inequality constraints is")
print(prob.constraints[0].dual_value)
