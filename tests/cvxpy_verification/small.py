import cvxpy as cp
import numpy as np

# Generate a random non-trivial linear program.

# Define and solve the CVXPY problem.
x = cp.Variable(2, nonneg=True)
prob = cp.Problem(cp.Minimize(-cp.sum(x)),
                  [cp.norm(cp.hstack([np.sqrt(6.) * x[0], x[1]])) <= 5, x >= 1])
prob.solve(solver=cp.ECOS, verbose=True)

# Print result.
print("\nThe optimal value is", prob.value)
print("A solution x is")
print(x.value)
print(np.sqrt((5**2 - 2) / 2))
# print("A dual solution is")
# print(prob.constraints[0].dual_value)
