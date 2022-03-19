import numpy as np

A = np.array([[0,1],
	     [-2,-3]])

x = np.array([[1],[1]])

def step(A, x):
	tmp = A@x
	tmp_norm = norm(tmp)
	return tmp/tmp_norm

def norm(x):
	return np.sqrt(np.sum(x ** 2))

step(A, x)

