import sys
sys.path.insert(0, '../dwavebinarycsp')
sys.path.insert(0, '../penaltymodel/penaltymodel_core')

import dwavebinarycsp
csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)
def any(a,b,c,d,e):
    return a or b or c or d or e
csp.add_constraint(any, ['a','b','c','d','e'])
bqm = dwavebinarycsp.stitch(csp)
