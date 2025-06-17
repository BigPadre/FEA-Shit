import numpy as np





'''
#part 1
#4 nodes: [0 1 2 3]
#3 elments: [0 1 2]
#element types = ['truss', 'truss', 'spring']
#n_DOF_Local_Truss = 2  # Number of local DOF for truss elements
#n_DOF_Local_Spring = 2  # Number of local DOF for spring elements
#n_DOF_global = 2  # Number of global DOF for each node
#theta_list = [120, 60, -90] deg (ordered by ascending node number)
#L_list = [5,5] m #length of each element 
#E_list = [200e9, 200e9] Pa #Young's modulus for each element
#Diameter_list = [0.025, 0.025] m #diameter of each element
#k_list=[400e3] N/m
#External_Loads = [0,-100e3] N global DOF ref list:[x,y]
#con_mat = [[0, 1], [0, 2], [0, 3]] 
#DOF_List_Global = [0, 1, 2, 3, 4, 5,6,7]  # Global DOF list for the nodes (2 GDOF*4 nodes)
#active_DOF= [0,1]#Node 0 global x and y DOFS
#fixed_DOF = [2:7]  
#initialize con mat
con_mat = np.array([[0, 1], [0, 2], [0, 3]])  # Connectivity matrix for the elements
#initilaize element geom and material properties
L_list = np.array([5, 5])  # Length of each element in meters
E_list = np.array([200e9, 200e9])  # Young's modulus for each element in Pascals
Diameter_list = np.array([0.025, 0.025])  # Diameter of each element in meters
k_list = np.array([400e3])  # Spring constant for the spring element in N/m
theta_list = np.array([120, 60, -90])*np.pi/180  # Angle of each element in degrees
A_list= np.pi * (Diameter_list / 2) ** 2  # Cross-sectional area for each element in m^2
#establish global and local DOFs
N_DOF_Local_Truss = 2  # Number of local DOF for truss elements
N_DOF_Local_Spring = 2  # Number of local DOF for spring elements
N_DOF_Global = 2  # Number of global DOF for each node
N_DOF_G_PER_ELEMENT = 4  # Number of global DOF per element (2 nodes * 2 DOF per node)
Tot_DOF= N_DOF_Global * (np.max(con_mat) + 1)  # Total degrees of freedom based on the highest node index
#fetch element types and their local stifness mats in terms of local DOFs (stifness factors calcd internally)
el_type_list=['truss','truss','spring']#ordered same as con mat

active_DOF = np.array([0, 1])  # Active DOF (global DOF indices)
fixed_DOF = np.array([2, 3, 4, 5, 6, 7])  # Fixed DOF (global DOF indices)
external_loads= np.array([0, -100e3, 0, 0, 0, 0, 0, 0])  # External loads at each node in global DOF format
'''
#part 3
#need beam shit bruh
#2D 3DOF/node 2 node els
#Global_DOF=3
#Total Nodes:4
#initialize con mat
con_mat = np.array([[0, 1],[1,2],[2,3]])  # Connectivity matrix for the elements


#initilaize element geom and material properties
L_list = np.array([240,240,240*np.sin(60*np.pi/180)])  # Length of each element in inches


E_list = np.array([29e6,29e6,29e6])  # Young's modulus for each element in Psi


r_o=3
r_i=2.5
I_beam=np.pi*((r_o**4)-(r_i**4))/4
I_list = np.array([I_beam,I_beam,I_beam])  # MOI of each element in in^4


k_list = np.array([None])  # Spring constant for the spring element in N/m

theta_list = np.array([90,0,120])*np.pi/180  # Angle of each element in degrees


A_beam=np.pi*((r_o**2)-(r_i**2)) #in^2


A_list= [A_beam,A_beam,A_beam]  # Cross-sectional area for each element in m^2
#establish global and local DOFs


N_DOF_Global = 3  # Number of global DOF for each node
N_DOF_G_PER_ELEMENT = 6  # Number of global DOF per element (2 nodes * 2 DOF per node)
Tot_DOF= N_DOF_Global * (np.max(con_mat) + 1)  # Total degrees of freedom based on the highest node index
#fetch element types and their local stifness mats in terms of local DOFs (stifness factors calcd internally)
el_type_list=['beam','beam','beam']#ordered same as con mat
DOF_list_Global = np.arange(Tot_DOF)  # Global DOF list for the nodes (3 GDOF*4 nodes)
fixed_DOF = np.array([0,1 ,9,10,11])  # Fixed DOF (global DOF indices x y t0,y1)

active_DOF = np.setdiff1d(DOF_list_Global, fixed_DOF)  # Active DOF (global DOF indices)

external_loads= np.array([0, 0, 0, 0, -2e3, 0])  # External loads at each node in global DOF format
#spring stiffness matrix
def spring_stiffness_matrix(k,N_DOF_G_PER_ELEMENT):
    if N_DOF_G_PER_ELEMENT==2:
        local_spring_stiffness_N_GDOF = np.array([[k, -k], [-k, k]])
    if N_DOF_G_PER_ELEMENT==4:
        local_spring_stiffness_N_GDOF = np.zeros((4, 4))
        local_spring_stiffness_N_GDOF[0, 0] = k
        local_spring_stiffness_N_GDOF[0, 2] = -k
        local_spring_stiffness_N_GDOF[2, 0] = -k
        local_spring_stiffness_N_GDOF[2, 2] = k
    return local_spring_stiffness_N_GDOF
#truss stiffness matrix
def truss_stiffness_matrix(L, E, A,N_DOF_G_PER_ELEMENT):
    k=(E * A / L)
    if N_DOF_G_PER_ELEMENT==2:
        local_truss_stiffness_N_GDOF = np.array([[k, -k], [-k, k]])
    if N_DOF_G_PER_ELEMENT==4:
        local_truss_stiffness_N_GDOF = np.zeros((4, 4))
        local_truss_stiffness_N_GDOF[0, 0] = k
        local_truss_stiffness_N_GDOF[0, 2] = -k
        local_truss_stiffness_N_GDOF[2, 0] = -k
        local_truss_stiffness_N_GDOF[2, 2] = k
    return local_truss_stiffness_N_GDOF
def beam_stiffness_matrix(L, E, I,A,N_DOF_G_PER_ELEMENT):
    # Local stiffness matrix for a 2D beam element with 6 DOF (2 translations and 1 rotation per node)
    c1=A*E/L
    c2=E*I/L**3
    if N_DOF_G_PER_ELEMENT==6:
        k_local_beam = np.array([
        [c1, 0, 0, -c1, 0, 0],
        [0, 12*c2, 6*c2*L, 0, -12*c2, 6*c2*L],
        [0, 6*c2*L, 4*c2*L**2, 0, -6*c2*L, 2*c2*L**2],
        [-c1, 0, 0, c1, 0, 0],
        [0, -12*c2, -6*c2*L, 0, 12*c2, -6*c2*L],
        [0, 6*c2*L, 2*c2*L**2, 0, -6*c2*L, 4*c2*L**2]
    ])
    elif N_DOF_G_PER_ELEMENT==4:
        k_local_beam=np.array([[12*c2, 6*c2*L,-12*c2, 6*c2],
                               [6*c2*L, 4*c2*L**2,-6*c2*L, 2*c2*L**2],
                               [-12*c2, -6*c2*L,12*c2, -6*c2],
                               [6*c2*L, 2*c2*L**2,-6*c2*L, 4*c2*L**2]])
    return k_local_beam

#convert local stiffness matrices to global stiffness matrices peices thru 4x4 2D transformation matrices (2 node elements, 2 Global DOFs per node)
def pick_T_mat (N_DOF_G_PER_ELEMENT,theta):
    if N_DOF_G_PER_ELEMENT == 4:
        T = np.array([[np.cos(theta), np.sin(theta), 0, 0],
                      [-np.sin(theta), np.cos(theta), 0, 0],
                      [0, 0, np.cos(theta), np.sin(theta)],
                      [0, 0, -np.sin(theta), np.cos(theta)]])
    elif N_DOF_G_PER_ELEMENT == 2:
        T = np.array([[np.cos(theta), np.sin(theta)],
                      [-np.sin(theta), np.cos(theta)]])
    elif N_DOF_G_PER_ELEMENT == 6:
        T = np.array([[np.cos(theta), np.sin(theta), 0, 0, 0, 0],
                      [-np.sin(theta), np.cos(theta), 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, np.cos(theta), np.sin(theta), 0],
                      [0, 0, 0, -np.sin(theta), np.cos(theta), 0],
                      [0, 0, 0, 0, 0, 1]])
    return T
#transform locals into global pieces
def Global_stiff_piece_generator(Local_stiffness_mat,theta,N_DOF_G_PER_ELEMENT):
    # Transform the local stiffness matrix to global coordinates
    T = pick_T_mat(N_DOF_G_PER_ELEMENT, theta)  # Transformation matrix
    K_global_piece = T.T @ Local_stiffness_mat @ T  # Transform the local stiffness matrix
    return K_global_piece
#assemble global stiffness matrix from global stiffness matrix pieces
def assemble_global_stiffness_matrix(k_global_piece_list, con_mat, Tot_DOF, N_DOF_Global):#might need to change N_DOF_element to a list of DOF per element
    
    Global_K = np.zeros((Tot_DOF, Tot_DOF))  # Initialize global stiffness matrix

    for i, k_global_piece in enumerate(k_global_piece_list):
        current_nodes = con_mat[i]  # e.g., [0, 1]
        # For 2 nodes, 2 DOF each: [0,1,2,3] for nodes 0 and 1
        dof_indices = []
        for node in current_nodes:
            dof_indices.extend([node * N_DOF_Global + d for d in range(N_DOF_Global)])
        Global_K[np.ix_(dof_indices, dof_indices)] += k_global_piece
    return Global_K
#fetch external loads and convert them to global DOF list format [[fx,fy],...len=N_nodes]


##assemble global load vector from distributed nodal loads and point loads
#indicate boundary conditions to reduce global stiffness matrix and global load vector

#solve for displacements using global stiffness matrix and global load vector
def solve_system(Global_K, Load_External_list, active_DOF):
    K_reduced = Global_K[np.ix_(active_DOF, active_DOF)]
    F_reduced = Load_External_list[active_DOF]
    displacements = np.linalg.solve(K_reduced, F_reduced)
    Tot_DOF = Global_K.shape[0]
    full_displacements = np.zeros(Tot_DOF)
    full_displacements[active_DOF] = displacements
    return full_displacements


def local_stiff_list(el_type_list):
    Spring_count=0
    Truss_count=0  
    beam_count=0
    Local_stiffness_mat_list = [None] * len(el_type_list)
    for i, el_type in enumerate(el_type_list):

        if el_type == 'spring':
            Local_stiffness_mat_list[i] = spring_stiffness_matrix(k_list[Spring_count],N_DOF_G_PER_ELEMENT) 
            Spring_count+=1

        elif el_type=='truss':
            Local_stiffness_mat_list[i] =truss_stiffness_matrix(L_list[Truss_count], E_list[Truss_count], A_list[Truss_count],N_DOF_G_PER_ELEMENT)
            Truss_count+=1
        elif el_type=='beam':
            Local_stiffness_mat_list[i] =beam_stiffness_matrix(L_list[beam_count], E_list[beam_count],I_list[beam_count], A_list[beam_count],N_DOF_G_PER_ELEMENT)
            beam_count+=1
    return Local_stiffness_mat_list
    
def Make_Global_piece_list(el_type_list,Local_stiffness_mat_list,theta_list,N_DOF_G_PER_ELEMENT):
    k_global_piece_list=[None] * len(el_type_list)
    for i in range(len(el_type_list)):
        k_global_piece_list[i] = Global_stiff_piece_generator(Local_stiffness_mat_list[i], theta_list[i],N_DOF_G_PER_ELEMENT)
    return k_global_piece_list
Local_stiffness_mat_list=local_stiff_list(el_type_list)
k_global_piece_list=Make_Global_piece_list(el_type_list,Local_stiffness_mat_list,theta_list,N_DOF_G_PER_ELEMENT)

Global_K=assemble_global_stiffness_matrix(k_global_piece_list, con_mat, Tot_DOF, N_DOF_G_PER_ELEMENT)
Deformation_list =solve_system(Global_K, external_loads, active_DOF)

print("Deformation Vector [x0 y0 t0 x1 y1 t1 ] m,m,rad:")
print(Deformation_list)