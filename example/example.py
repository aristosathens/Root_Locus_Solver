from root_locus_solver import root_locus
from pprint import pprint

poles = [ -4 + 2j, -4 - 2j, 1 ]
zeros = [ -1 ]

values = root_locus(poles = poles, zeros = zeros)
pprint(values)