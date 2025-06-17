import numpy as np


#constant strain triangle displacement field and stresses calculator
# Total degrees of freedom for a constant strain triangle (2 nodes, 3 DOF each)
DOF_list= np.array([0, 1, 2, 3, 4, 5])  # Global degrees of freedom for a constant strain triangle
active_DOF=np.array([4,5])
fixed_DOF = np.setdiff1d(DOF_list, active_DOF)

def solve_system(Global_K, Load_External_list, active_DOF):
    K_reduced = Global_K[np.ix_(active_DOF, active_DOF)]
    F_reduced = Load_External_list[active_DOF]
    displacements = np.linalg.solve(K_reduced, F_reduced)
    Tot_DOF = Global_K.shape[0]
    full_displacements = np.zeros(Tot_DOF)
    full_displacements[active_DOF] = displacements
    return full_displacements
#const strain tri for wide ass part: 6x6 stiffness mat
def const_strain_tri_disp_field_and_stresses(x_list, y_list,ratio,E,active_DOF=None,Force=None,thin=False,  u_node_list=None, v_node_list=None,t=None):
    # Calculate twice the area of the triangle
    two_A = np.dot(x_list, [y_list[1]-y_list[2], y_list[2]-y_list[0], y_list[0]-y_list[1]])
    # Shape function derivatives
    b = np.array([y_list[1] - y_list[2], y_list[2] - y_list[0], y_list[0] - y_list[1]]) / two_A# betas
    c = np.array([x_list[2] - x_list[1], x_list[0] - x_list[2], x_list[1] - x_list[0]]) / two_A#gammas
    strains=np.zeros(3)
    disp_vector=np.zeros(6)
    #build B mat out of gammas and betas
    B=np.array([[b[0], 0, b[1], 0, b[2], 0],
              [0, c[0], 0, c[1], 0, c[2]],
              [c[0], b[0], c[1], b[1], c[2], b[2]]])
    
        # Displacement field at each node (displ applied or nah)
    if u_node_list or v_node_list!=None:
        u_disp_field = np.array(u_node_list)
        v_disp_field = np.array(v_node_list)
        disp_vector= [u_disp_field, v_disp_field]
        strains= B@disp_vector #ex, ey, exy
        stresses = D @ strains #sx, sy, sxy
    else:
        u_disp_field = np.zeros(3)
        v_disp_field = np.zeros(3)
    #build D mat out of E and ratio for stiff mat, or just k mat depending on thickness conidtion 
    if thin== True:
        if t!=None:
            print("Thin or not?")
            return
    #thin case
        D= np.array([[E/(1-ratio**2), E*ratio/(1-ratio**2), 0],
              [E*ratio/(1-ratio**2), E/(1-ratio**2), 0],
              [0, 0, E/(2*(1+ratio))]])
        
    else:    
        D=np.array([[1-ratio,ratio,0],[ratio,1-ratio,0],[0,0,.5-ratio]])*(E/((1+ratio)*(1-2*ratio)))
    
    k_local=B.T@D@B*t*.5*two_A  # totally didn't forget to multiply by A to match the final k  in the book
    # If Force is provided, calculate the displacement field based on the force
    if Force is not None:
        # Assuming Force is a vector of nodal forces, calculate displacements
        if len(Force) != 6:
            raise ValueError("Force vector must have 6 components for a 2D constant strain triangle.")
    disp_vector = solve_system(k_local, Force, active_DOF)
    # Calculate strains using the B matrix and the displacement vector
    strains = B @ disp_vector  # ex, ey, exy

    # Calculate stresses using the B matrix and D matrix, strains
    stresses = D @ strains  # sx, sy, sxy

    return disp_vector, strains,stresses,k_local


disp_field, strains,stresses,k_local = const_strain_tri_disp_field_and_stresses(np.array([0,1,0.5]), np.array([0,0,1.25]), 0.25,30e6,active_DOF,np.array([0,0,0,0,1.5e6,0]),thin=False,u_node_list=None, v_node_list=None,t=24)

print("\n k_local:")
print(k_local)
print("disp field [x0,y0,x1,y1,x2,y2]:")
print(disp_field)
print("\n strains [ex, ey, exy]:")
print(strains)
print("\n stresses [sx, sy, sxy]:")
print(stresses)