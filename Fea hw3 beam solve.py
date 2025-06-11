import numpy as np

# Problem parameters
I=16.014#in^4
L=72#in
E=29e6#psi
F=1000#lbf
N_e = 1  # Number of elements
N_DOF = 2  # Degrees of freedom per node
N_Nodes = 2  # Number of nodes
Tot_DOF = N_DOF * N_Nodes  # Total degrees of freedom
node_list = np.array(list(range(0, N_Nodes)))  # 0 index convention
full_displacements = np.zeros(Tot_DOF)  # Initialize full displacement vector
current_nodes = np.zeros([1, N_DOF])  # Node list for each element
'''
defl_max=-F*L**3/(3*E*I) #max deflection at free end
print(f"Theoretical deflection at free end: {defl_max:.4f} in")
'''
def beam_stiffness_mat(L,E,I):
    return np.array([[12,6*L,-12,6*L],[6*L,4*L**2,-6*L,2*L**2],[-12,-6*L,12,-6*L],[6*L,2*L**2,-6*L,4*L**2]])*(E*I/L**3)
def Global_T_Mat(theta):
    return np.array([
          [ np.cos(theta),  np.sin(theta), 0, 0],
    [-np.sin(theta),  np.cos(theta), 0, 0],
    [ 0,  0, np.cos(theta), np.sin(theta)],
    [ 0,  0,-np.sin(theta), np.cos(theta)]])
def assemble_global_stiffness_matrix(k_global_piece_list, con_mat, N_DOF):
    Tot_DOF = 2 * (np.max(con_mat) + 1)  # Total degrees of freedom based on the highest node index
    Global_K = np.zeros((Tot_DOF, Tot_DOF))  # Initialize global stiffness matrix

    for i, k_global_piece in enumerate(k_global_piece_list):#i assume that enumerate gives the indexed matrix of the element in the piece list?
        current_nodes = con_mat[i]  # Get the nodes for the current element
        dof_indices = [current_nodes[0] * N_DOF, current_nodes[0] * N_DOF + 1,
                       current_nodes[1] * N_DOF, current_nodes[1] * N_DOF + 1]
        Global_K[np.ix_(dof_indices, dof_indices)] += k_global_piece

    return Global_K
def assemble_applied_load_vector(num_elements, con_mat, dist_nodal_loads=None,point_loads=None):
    # Initialize the global load vector
    Tot_DOF = 2 * (np.max(con_mat) + 1)  # Total degrees of freedom based on the highest node index
    global_load_vector = np.zeros(Tot_DOF)

    # Add distributed loads
    if dist_nodal_loads is not None:
        for i in range(num_elements):
            current_nodes = con_mat[i]  # Get the nodes for the current element
            dof_indices = [current_nodes[0] * 2, current_nodes[0] * 2 + 1, current_nodes[1] * 2, current_nodes[1] * 2 + 1]
            global_load_vector[dof_indices] += dist_nodal_loads[i]

    # Add point loads
    if point_loads is not None:
        for node_index, load_value in point_loads.items():
            dof_index = node_index * 2  # Assuming point loads are applied in the x-direction
            global_load_vector[dof_index] += load_value

    return global_load_vector
def construct_and_transform_local_stiffness_matrices(num_elements, type_el_list, L_list, E_list, I_list, k_list,theta_list):
    
    k_global_piece_list = []  # List to store global stiffness matrix pieces
    for i in range(num_elements):
        current_element = type_el_list[i]  # Get the type of the current element
        if current_element == 'beam':
            L = L_list[i]
            E = E_list[i]
            I = I_list[i]
            K_local_mat = beam_stiffness_mat(L, E, I)  # or spring_stiffness_mat(k)
        elif current_element == 'spring':
            k = k_list[i]
            K_local_mat = spring_stiffness_mat(k)
        else:
            raise ValueError(f"Unknown element type: {current_element}")
        
        # Transform the local stiffness matrix to global coordinates
        theta = theta_list[i]  # Get the angle for the current element
        T = Global_T_Mat(theta)  # Transformation matrix
        K_global_piece = T.T @ K_local_mat @ T  # Transform the local stiffness matrix
        k_global_piece_list.append(K_global_piece)  # Append the global stiffness pieces  to the list


    return k_global_piece_list#list of global stiffness matrices ordered by element number
def calculate_stresses_in_elements(k_global_piece_list, con_mat, full_displacements, N_DOF):#stress list based on global sitfnesses and displacements of nodes and els
#Calculate Stresses in each element
    stress_list = np.zeros(len(k_global_piece_list))  # Initialize stress list for each element
    for i in range(len(k_global_piece_list)):
        k_local = k_global_piece_list[i]
        current_nodes = con_mat[i]  # Get the nodes for the current element
        dof_indices = [current_nodes[0] * N_DOF, current_nodes[0] * N_DOF + 1,
                       current_nodes[1] * N_DOF, current_nodes[1] * N_DOF + 1]
        # Calculate stress in the element with forces along element axis
        element_displacements = full_displacements[dof_indices]
        element_stress = k_local @ element_displacements  # Calculate stress using local stiffness matrix
        stress_list[i] = element_stress[0]  # Assuming the first component is the stress in the element
     
    return stress_list
def solve_system(Global_K, Load_External_list, active_DOF, k_global_piece_list, con_mat, N_DOF):#will NDOF be used correctly?
    # Reduce the system for active DOF
    K_reduced = Global_K[np.ix_(active_DOF, active_DOF)]
    F_reduced = Load_External_list[active_DOF]

    # Solve for displacements at active degrees of freedom
    displacements = np.linalg.solve(K_reduced, F_reduced)
    # Initialize the full displacement vector
    Tot_DOF = Global_K.shape[0] 
    full_displacements = np.zeros(Tot_DOF)
    # Reconstruct the full displacement vector
    full_displacements[active_DOF] = displacements
    stresses= calculate_stresses_in_elements(k_global_piece_list, con_mat, full_displacements, N_DOF)
    return displacements,stresses
con_mat= np.array([[0, 1]])  # Connectivity matrix for the truss (element connects node 0 and node 1)
type_el_list = ['beam']  # List of element types
L_list = [L]  # Length of each element
E_list = [E]  # Young's modulus for each element
I_list = [I]  # Moment of inertia for each element
num_elements = len(con_mat)  # Number of elements
k_list = [None] * num_elements  # Placeholder for spring stiffness, not used in this case
theta_list = [0]  # Angles for each element in radians (0 for horizontal beam)
# Construct and transform local stiffness matrices
k_global_piece_list = construct_and_transform_local_stiffness_matrices(num_elements, type_el_list, L_list, E_list, I_list, k_list, theta_list)
# Assemble global stiffness matrix
Global_K = assemble_global_stiffness_matrix(k_global_piece_list, con_mat, N_DOF=2)
# Assemble applied load vector
# Create the point_loads dictionary for node 1, y-direction (DOF 3)
point_loads = {1: -F}  # This assumes your function applies the load to the correct DOF

Load_External_Vector = assemble_applied_load_vector(
    num_elements, con_mat, dist_nodal_loads=None, point_loads=point_loads
)

active_DOF = [2,3]  
fixed_DOF = [0, 1]
# Solve the system
displacements, stresses = solve_system(Global_K, Load_External_Vector, active_DOF, k_global_piece_list, con_mat, N_DOF=2)
# Print results
print(f"Displacement at node 1: {displacements[0]:.4f} in")
print(f"Stress in element 0: {stresses[0]:.4f} psi")
reactions = Global_K[np.ix_(fixed_DOF, range(Tot_DOF))] @ full_displacements - Load_External_Vector[fixed_DOF]
