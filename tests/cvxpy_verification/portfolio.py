import cvxpy as cp
import numpy as np

gamma = 0.5
F = np.array([0.604897, 0.0452059,
              0.329554, 0.257742,
              0.536459, 0.270431,
              0.444451, 0.0268018,
              0.10794, 0.904459]).reshape(5, 2)
Sigma = np.array([1.20033, 0.210998, 0.336728, 0.270059, 0.106179, 0.210998, 0.44646, 0.246494, 0.153379, 0.268689, 0.336728, 0.246494, 0.795515,
                  0.245678, 0.302499, 0.270059, 0.153379, 0.245678, 0.91505, 0.0722151, 0.106179, 0.268689, 0.302499, 0.0722151, 1.04364]).reshape(5, 5)
mu1 = np.array([0.680375,
                0.211234,
                0.566198,
                0.59688,
                0.823295])
mu2 = np.array([0.967399,
                0.514226,
                0.725537,
                0.608354,
                0.686642])

# Define and solve the CVXPY problem.
x = cp.Variable(5, nonneg=True)
prob = cp.Problem(cp.Minimize(cp.quad_form(x, gamma*Sigma) - mu1 * x),
                  [x >= 0, cp.sum(x) == 1])
prob.solve(solver=cp.OSQP, verbose=True)
# [0.24424712, 0., 0.01413456, 0.25067381, 0.4909445 ]

x = cp.Variable(5, nonneg=True)
prob = cp.Problem(cp.Minimize(cp.quad_form(x, gamma*Sigma) - mu2 * x),
                  [x >= 0, cp.sum(x) == 1])
prob.solve(solver=cp.OSQP, verbose=True)
# [4.38579051e-01, 3.04375987e-23, 2.00025310e-01, 1.17002001e-01, 2.44393639e-01]

# Print result.
print("\nThe optimal value is", prob.value)
print("A solution x is")
print(x.value)
# print("A dual solution is")
# print(prob.constraints[0].dual_value)
