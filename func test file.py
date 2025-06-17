
'''
def convert_dist_loads_to_nodal_loads(elements_containing_dist_loads, load_type_list, con_mat, Tot_DOF):#global vector of distributed loads converted to nodal loads
    # Convert distributed loads to nodal loads
    dist_nodal_loads_list = np.zeros(Tot_DOF)  # Initialize the distributed nodal loads vector
    for i, element in enumerate(elements_containing_dist_loads):
        current_nodes = con_mat[element]
        for j in range(len(current_nodes)):
            local_dof_indices = [current_nodes[j] * N_DOF_Nodes + k for k in range(N_DOF_Nodes)]
            if load_type_list[i] >0:
                dist_nodal_loads_list[local_dof_indices]
      # Distribute load evenly among nodes
    return dist_nodal_loads_list
    '''
#func to convert solved global defs to local defs to calculate stresses in elements
#input: def vector, el list, con mat, local stifness mats, Alist
def calc_stresses_per_el(el_type_list,N_DOF_G_PER_ELEMENT,A_list,loc_def_vector):
    #index el type list 
    
     
    #fetch local stiff mat in global DOF form using external if
        
    #fetch node indexes per el (index con mat from type list number)
    #fetch defs acting on node from global perspective
    #transform defs into local coords
    local_def_list
    #multiply the transformed defs per element by prev local k mat for force developed
    Force_dev_list=local_stiffness_list(el_type_list) @ loc_def
    Stress_list=Force_dev_list/A_list
    #Filter out if spring element cause ts don't have a crosssection 
    #divide by A list indexed by el type index

def Global_to_local_Transform (theta,global_def,):
