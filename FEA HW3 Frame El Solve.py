import numpy as np
# Problem parameters
N_e = 2  # Number of elements
N_DOF = 2  # Degrees of freedom per node
N_Nodes = 3  # Number of nodes
Tot_DOF = N_DOF * N_Nodes  # Total degrees of freedom
node_list=np.array(list(range(0,N_Nodes)))#0 index convention
current_nodes = np.zeros([1,N_DOF])#node list for 
k_loc_list = []  # List to store local stiffness matrices

# Boundary conditions
active_DOF = [2]  # DOFs for node 1 (x)
fixed_DOF = [i for i in range(Tot_DOF) if i not in active_DOF]
# Connectivity matrix and local stiffness values


con_mat = np.array([[node_list[0],node_list[1]], [node_list[1],node_list[2]]])  # Connectivity matrix 0 index convention

A_list = np.array([0.03**2,0.01**2])*(np.pi/4)  # Cross-sectional area (m^2)
E_list = [200e9,100e9]  # Young's modulus (Pa)
L_list=[1,0.33]#m
# Calculate local stiffness values for each element

for i in range(N_e):

    k_loc_list.append(A_list[i] * E_list[i] / L_list[i])  # Stiffness for each element




# Initialize global stiffness matrix
Global_K = np.zeros((Tot_DOF, Tot_DOF))

# Function to evaluate the global transformation matrix
def Global_T_Mat(theta):
    return np.array([
        [np.cos(theta)**2, np.sin(theta)*np.cos(theta), -np.cos(theta)**2, -np.sin(theta)*np.cos(theta)],
        [np.sin(theta)*np.cos(theta), np.sin(theta)**2, -np.sin(theta)*np.cos(theta), -np.sin(theta)**2],
        [-np.cos(theta)**2, -np.sin(theta)*np.cos(theta), np.cos(theta)**2, np.sin(theta)*np.cos(theta)],
        [-np.sin(theta)*np.cos(theta), -np.sin(theta)**2, np.sin(theta)*np.cos(theta), np.sin(theta)**2]
    ])

# Function to assemble truss stiffness matrices
def Stiff_Local_to_Global(k, theta):
    return k * Global_T_Mat(theta)

# Populate global stiffness matrix
thetas = [np.radians(0), np.radians(0)]  # Angles for each element in radians


for i in range(N_e):
    
    theta = thetas[i]  # Get the angle for the current element
    current_nodes = con_mat[i][0], con_mat[i][1]#node list per el from con mat
    
    dof_indices = [current_nodes[0]*N_DOF, current_nodes[0]*N_DOF+1, current_nodes[1]*N_DOF, current_nodes[1]*N_DOF+1] # Expand indices to global DOF indicies (per element)
    # Add the element stiffness matrix to the global stiffness matrix
    Global_K[np.ix_(dof_indices, dof_indices)] += Stiff_Local_to_Global(k_loc_list[i], thetas[i])




# Applying force
Force_External = np.zeros(Tot_DOF)

# Define the magnitude and angle of the force applied at node 3
F_magnitude = 55e3  # Force magnitude in N
F_angle = np.radians(0)  # Force angle in degrees (convert to radians)

# Decompose the force into x and y components
F_x = F_magnitude * np.cos(F_angle)
F_y = F_magnitude * np.sin(F_angle)

# Apply the force components to the DOFs of node 1
Force_External[2] = F_x  # x-component of the force at DOF 2


# Reduce the system for active DOF
K_reduced = Global_K[np.ix_(active_DOF, active_DOF)]
F_reduced = Force_External[active_DOF]

# Solve for displacements at active degrees of freedom
disp_reduced = np.linalg.solve(K_reduced, F_reduced)

# Reconstruct the full displacement vector (Units of meters)
disp_list = np.zeros(Tot_DOF)
disp_list[active_DOF] = disp_reduced

stress_list = np.zeros(N_e)  # Initialize stress list for each element
# Calculate stress in each element
for i in range(N_e):
    # Get the local stiffness for the current element
    k_local = k_loc_list[i]
    
    # Get the DOF indices for the current element
    current_nodes = con_mat[i]
    dof_indices = [current_nodes[0]*N_DOF, current_nodes[0]*N_DOF+1, current_nodes[1]*N_DOF, current_nodes[1]*N_DOF+1]
    
    # Calculate the force in the element using dot global stiffness matrix and displacements
    force_element = np.dot(Global_K[np.ix_(dof_indices, dof_indices)], disp_list[dof_indices])
    
    # Calculate stress in the element (stress = force / area)
    stress_list[i] = force_element[0] / A_list[i]  # Assuming uniform stress across the element
for i in range(N_e):
    current_nodes = con_mat[i]
    dof_indices = [
        current_nodes[0]*N_DOF, current_nodes[0]*N_DOF+1,
        current_nodes[1]*N_DOF, current_nodes[1]*N_DOF+1
    ]
    # Get global displacements for the element
    u_global = disp_list[dof_indices]
    # Get element angle
    theta = thetas[i]
    # Transformation for local axis
    c = np.cos(theta)
    s = np.sin(theta)
    # Local displacement difference (axial)
    u1_local = c * u_global[0] + s * u_global[1]
    u2_local = c * u_global[2] + s * u_global[3]
    axial_disp = u2_local - u1_local
    # Axial force and stress
    axial_force = k_loc_list[i] * axial_disp
    stress_list[i] = axial_force / A_list[i]
Reaction_forces = np.zeros(Tot_DOF)  # Initialize reaction forces vector
# Calculate reaction forces at fixed DOFs
for i in fixed_DOF:
    Reaction_forces[i] = np.dot(Global_K[i, :], disp_list)  # Reaction force at fixed DOF
# Output results
print("Global Stiffness Matrix (N/m):")
print(Global_K)
print("\nDisplacement Vector (m) [x0 y0 x1 y1 x2 y2]:")
print(disp_list)
print("\nReaction Forces at Fixed DOFs (N)[x0 y0 x1 y1 x2 y2]:")
print(Reaction_forces)
print("\nStress in each element (Pa):")
print(stress_list)