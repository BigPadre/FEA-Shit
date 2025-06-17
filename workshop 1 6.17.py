import numpy as np


def const_strain_tri_disp_field_and_stresses(x_list, y_list, u_node_list, v_node_list,ratio,E,thin):
    # Calculate twice the area of the triangle
    two_A = np.dot(x_list, [y_list[1]-y_list[2], y_list[2]-y_list[0], y_list[0]-y_list[1]])
    # Shape function derivatives
    b = np.array([y_list[1] - y_list[2], y_list[2] - y_list[0], y_list[0] - y_list[1]]) / two_A# betas
    c = np.array([x_list[2] - x_list[1], x_list[0] - x_list[2], x_list[1] - x_list[0]]) / two_A#gammas
    #build B mat out of gammas and betas
    B=np.array([[b[0], 0, b[1], 0, b[2], 0],
              [0, c[0], 0, c[1], 0, c[2]],
              [c[0], b[0], c[1], b[1], c[2], b[2]]])
    #build D mat out of E and ratio 
    if thin== True:
    #thin
        D= np.array([[E/(1-ratio**2), E*ratio/(1-ratio**2), 0],
              [E*ratio/(1-ratio**2), E/(1-ratio**2), 0],
              [0, 0, E/(2*(1+ratio))]])
    else:
        D=np.array([[1-ratio,ratio,0],[ratio,1-ratio,0],[0,0,.5-ratio]])*(E/((1+ratio)*(1-2*ratio)))
    # Strain components
    ex = np.sum(b * u_node_list)
    ey = np.sum(c * v_node_list)
    exy = np.sum(b * v_node_list) + np.sum(c * u_node_list)

    # Displacement field at each node (list not array input)
    u_disp_field = np.array(u_node_list)
    v_disp_field = np.array(v_node_list)
    def_field_const_strain_tri = [u_disp_field, v_disp_field]
    strains= np.array([ex, ey, exy])

    # Calculate stresses using the B matrix and D matrix, strains
    stresses = D @ strains #sx, sy, sxy
    return def_field_const_strain_tri, strains,stresses


disp_field, strains,stresses = const_strain_tri_disp_field_and_stresses([0,2,0], [-1,0,1], [0,0.0012,0], [.0025,0,.0025],0.25,30e6,True)
print("disp field [x0,x1,x2],[y0,y1,y2]:")
print(disp_field)
print("\n strains [ex, ey, exy]:")
print(strains)
print("\n stresses [sx, sy, sxy]:")
print(stresses)