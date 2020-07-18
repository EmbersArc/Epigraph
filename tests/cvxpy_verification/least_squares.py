# Import packages.
import cvxpy as cp
import numpy as np

# Generate data.
m = 20
n = 15
np.random.seed(1)
A = np.random.randn(m, n) * (np.random.randn(m, n) > 0.)
b = np.random.randn(m)

print(repr(A.reshape(1, -1)))
print(repr(b))

# Define and solve the CVXPY problem.
x = cp.Variable(n)
cost = cp.sum_squares(A @ x - b)
prob = cp.Problem(cp.Minimize(cost))
prob.solve(solver=cp.OSQP, verbose=True)

# Print result.
print("\nThe optimal value is", prob.value)
print("The optimal x is")
print(repr(x.value))
print("The norm of the residual is ", cp.norm(A @ x - b, p=2).value)