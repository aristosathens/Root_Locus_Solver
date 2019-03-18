from root_locus_solver import root_locus, print_info

poles = [ -4 + 2j, -4 - 2j, 1 ]
zeros = [ -1 ]

values = root_locus(poles = poles, zeros = zeros)
print_info(values)

'''
    Output:
    
    ------ <module> ------
    values:
    {'asymptote_angles': array([-90.,  90.]),
    'centroid': (-3+0j),
    'departure_angles': [281.88865803962796, 78.11134196037204, 180.0],
    'poles': [(-4+2j), (-4-2j), 1],
    'real_axis_angles': [],
    'real_axis_points': array([], dtype=complex128),
    'zeros': [-1]}
'''