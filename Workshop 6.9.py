import numpy as np

A = 15.55e5 * np.array([72, -18, 18])
B = np.array([-18*15.55e5, (200000)+12*15.55e5 , 18*15.55e5])
C = 15.55e5 * np.array([18, -18, 36])


K = np.vstack([A, B, C])

F= np.array([0, -50e3, 0]) #external loads
# Solve for displacements
disp = np.linalg.solve(K, F)
print("Displacements: [theta2, v3, theta 3]", disp)