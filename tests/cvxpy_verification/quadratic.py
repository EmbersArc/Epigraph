import cvxpy as cp
import numpy as np

# Define and solve the CVXPY problem.
x = cp.Variable(3)
P = np.array([[1, 0.5, 0], [0.5, 1, 0], [0, 0, 2]])
c = np.array([1, 3, 2])

obj = cp.quad_form(x, P) + c @ x
prob = cp.Problem(cp.Minimize(obj),
                  [cp.sum(x) == 1, -1 <= x, x <= 1])
prob.solve(solver=cp.OSQP, verbose=True)

# Print result.
print("\nThe optimal value is", prob.value)
print("A solution x is")
print(x.value)
# print("A dual solution is")
# print(prob.constraints[0].dual_value)
