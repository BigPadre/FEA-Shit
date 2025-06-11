import numpy as np

#code for 6 DOF system, needs to be translated for a 4 DOF system with triangular loads

#2d stiff mat for beam el (x y  translations and rotation about z) 2nodes: 6x6 local
def Two_D_6DOF_Global_peice_Transformed_list(el_type_list,L_list, E_list, I_list,theta_list,A_list,k_list):
# Transform the local stiffness matrix to global coordinates
    k_global_piece_list = []  # List to store global stiffness matrix pieces
    for i in range(len(el_type_list)):
        if el_type_list[i]== 'beam':
            L = L_list[i]
            E = E_list[i]
            I = I_list[i]
            
            A=A_list[i]
            K_local_mat = Local_Two_D_6DOF_beam_stiffness_mat(L, E, I,A)  # Local stiffness matrix for a beam element
            
            K_global_piece = T.T @ K_local_mat @ T  # Transform the local stiffness matrix
        elif el_type_list[i] == 'spring':
            k=k_list[i]
            K_global_piece = np.array([[k, 0, 0, -k, 0, 0],
                                       [0, 0, 0, 0, 0, 0],
                                       [0, 0, 0, 0, 0, 0],
                                       [-k, 0, 0, k, 0, 0],
                                       [0, 0, 0, 0, 0, 0],
                                       [0, 0, 0, 0, 0, 0]])
        theta = theta_list[i]
        T = Two_D_6DOF_T_mat(theta)  # Transformation matrix
        k_global_piece_list.append(K_global_piece)  # Append the global stiffness pieces  to the list
    return k_global_piece_list

def Local_Two_D_6DOF_beam_stiffness_mat(L, E, I,A):
    # Local stiffness matrix for a 2D beam element with 6 DOF (2 translations and 1 rotation per node)
    c1=A*E/L
    c2=E*I/L**3
    k_local = np.array([[c1, 0, 0, -c1, 0, 0],
                          [0, 12*c2, 6*c2*L, 0, -12*c2, 6*c2],
                          [0, 6*c2*L, 4*c2*L**2, 0, -6*c2*L, 2*c2*L**2],
                          [-c1, 0, 0, c1, 0, 0],
                          [0, -12*c2, -6*c2*L, 0, 12*c2, -6*c2],
                          [0, 6*c2*L, 2*c2*L**2, 0, -6*c2*L, 4*c2*L**2]])
    return k_local
def Two_D_6DOF_T_mat(theta):
    return np.array([[ np.cos(theta),  np.sin(theta), 0, 0, 0, 0],
    [-np.sin(theta),  np.cos(theta), 0, 0, 0, 0],
    [ 0,  0, 1, 0,0 , 0],
    [0,0, 0, np.cos(theta), np.sin(theta), 0],
    [0,0, 0,-np.sin(theta), np.cos(theta), 0],
    [ 0,  0,0, 0, 0, 1]])
    

def assemble_applied_load_vector(num_elements, con_mat, dist_nodal_loads=None,point_loads=None,Tot_DOF=None):
    # Initialize the global load vector
    global_load_vector = np.zeros(Tot_DOF)
    if Tot_DOF is None:
        print("Dumbass")
    # Add distributed loads
    if dist_nodal_loads is not None:
        for i in range(num_elements):
            current_nodes = con_mat[i]  # Get the nodes for the current element
            dof_indices = [current_nodes[0] * 2, current_nodes[0] * 2 + 1, current_nodes[1] * 2, current_nodes[1] * 2 + 1]
            global_load_vector[dof_indices] += dist_nodal_loads[i]

    # Add point loads
    if point_loads is not None:
        global_load_vector += point_loads
    return global_load_vector

def assemble_global_stiffness_matrix(k_global_piece_list, con_mat, Tot_DOF, N_DOF_element):#might need to change N_DOF_element to a list of DOF per element
    
    Global_K = np.zeros((Tot_DOF, Tot_DOF))  # Initialize global stiffness matrix

    for i, k_global_piece in enumerate(k_global_piece_list):#i assume that enumerate gives the indexed matrix of the element in the piece list?
        current_nodes = con_mat[i]  # Get the nodes for the current element (hopefully assigns them automatically)
        for j in range(len(current_nodes)):
            dof_indices = [current_nodes[j] * N_DOF_element + k for k in range(N_DOF_element)]
        
        Global_K[np.ix_(dof_indices, dof_indices)] += k_global_piece
    return Global_K

def solve_system(Global_K, Load_External_list, active_DOF):
    K_reduced = Global_K[np.ix_(active_DOF, active_DOF)]
    F_reduced = Load_External_list[active_DOF]
    displacements = np.linalg.solve(K_reduced, F_reduced)
    Tot_DOF = Global_K.shape[0]
    full_displacements = np.zeros(Tot_DOF)
    full_displacements[active_DOF] = displacements
    return full_displacements


L_list=np.array([10,10,10])*12#inches
E_list=np.array([30e6, 30e6, 30e6])  # Young's modulus in psi   
I_list=np.array([200, 100, 200])  # Moment of inertia in in^4
A_list=np.array([10,10,10])  # Cross-sectional area in in^2
theta_list=np.array([np.pi/2, 0, -np.pi/2])  # Angles in radians
num_beam_el = len(L_list)  # Number of beam elements
con_mat = np.array([[0, 1], [1, 2], [2, 3]])  # Connectivity matrix for the beam elements
n_dof_beam_nodes=3
point_loads = np.array([0,0,0,10e3,0,0,0,0,5e3,0,0,0])  # Point load array [[fx fy m]...Tot_DOF for indexes]
dist_nodal_loads = None  # Distributed nodal loads array [f m...N_elements]
N_DOF_Nodes=n_dof_beam_nodes  # Degrees of freedom per node (3 for 2D beam: x, y, rotation)
Tot_DOF = N_DOF_Nodes * (np.max(con_mat) + 1)  # Total degrees of freedom based on the highest node index
fixed_DOF = np.array([0, 1, 2, 9, 10,11])  # Active degrees of freedom (indices of the full displacement vector)
active_DOF = np.setdiff1d(np.arange(Tot_DOF), fixed_DOF)  # Active DOF are those not fixed
Global_pieces= Two_D_6DOF_beam_stiffness_mat_Global_peice(num_beam_el, L_list, E_list, I_list, theta_list, A_list)  # Construct global stiffness matrix pieces
Global_K = assemble_global_stiffness_matrix(Global_pieces, con_mat, Tot_DOF, N_DOF_Nodes)
Load_External_list = assemble_applied_load_vector(num_beam_el, con_mat, dist_nodal_loads, point_loads, Tot_DOF)
full_displacements = solve_system(Global_K, Load_External_list, active_DOF)
reactions = Global_K[np.ix_(fixed_DOF, range(Tot_DOF))] @ full_displacements - Load_External_list[fixed_DOF]
# Print results
print("Displacement Vector:")
print(full_displacements)
print("\nReactions at fixed DOF [fx0 fy0 Mz0 fx3 fy3 Mz3]:")
print(reactions)