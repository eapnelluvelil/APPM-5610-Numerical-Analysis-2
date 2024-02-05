import numpy as np
import scipy.special

def odefun(t, y_0):
	if (t == 0):
		return np.array([0, -1/8])
	else:
		return np.array([y_0[1], (-3/t)*y_0[1] - y_0[0]])

def trap(odefun, t_0, t_f, h, y_0):
	y = y_0
	while (t_0 <= t_f):
		E_y = y + (h/2)*odefun(t_0, y)
		I_m = np.array([[1, (-h/2)], [(h/2), (1 + (3/2)*(h/(t_0 + h)))]])
		I_y = np.linalg.solve(I_m, E_y)

		y = E_y + (h/2)*odefun(t_0 + h, I_y)
		t_0 = t_0 + h
	return y

t_0 = 0
t_f = 3*np.arccos(-1)
h = (t_f - t_0)
y_0 = np.array([1/2, 0])
tol = 10**(-11)

num_soln = 0

max_rows = 21
A = np.zeros((max_rows, max_rows))

y = trap(odefun, t_0, t_f, h, y_0)
A[0, 0] = y[0]

for i in range(1, max_rows):
	h /= 2
	y = trap(odefun, t_0, t_f, h, y_0)
	A[i, 0] = y[0]

	for j in range(1, i+1):
		A[i, j] = ((4**j)*A[i, j-1] - A[i-1, j-1])/(4**j - 1)

	if (np.abs(A[i, i] - A[i-1, i-1]) < tol):
		num_soln = A[i, i]
		break

actual_soln = scipy.special.jv(1, t_f) 

print(np.abs(actual_soln - t_f*num_soln))
print(np.abs(actual_soln - t_f*A[max_rows-1, max_rows-1]))