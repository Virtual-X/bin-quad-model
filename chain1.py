import dwavebinarycsp
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx
import matplotlib.pyplot as plt

import dimod

def createNames(base, count):
    return [base+str(i)+'_' for i in range(count)]

def createToupleSet(nval):
    return { ((0,)*i) + (1,) + ((0,)*(nval - i - 1)) for i in range(nval) }

links = [(1,2), (2,0), (2,0), (0,1), (1,2), (0,1)]

nlinks = len(links)
linknames = createNames('link', nlinks)
linkbits = createToupleSet(nlinks)

nnodes = nlinks - 1
nodenames = createNames('nodes', nnodes)
nodebits = 2

csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

def not_all_2(a, b):
    return not (a and b)

for name in linknames:
    variables = [name+str(i) for i in range(nlinks)]
    csp.add_constraint(linkbits, variables)

for i in range(nlinks):
    variables = [name+str(i) for name in linknames]
    csp.add_constraint(linkbits, variables)

def create_if_a_then_bc_maps_to(v):
    v0 = v & 1
    v1 = v >> 1
    def if_a_then_bc_maps_to_v(a, b, c):
        return not a or (b == v0 and c == v1)
    return if_a_then_bc_maps_to_v

for k in range(nlinks):
    n = linknames[k]
    if k > 0:
        l = nodenames[k-1]
        for i in range(nlinks):
            vars = [n+str(i), l+'0', l+'1']
            csp.add_constraint(create_if_a_then_bc_maps_to(links[i][0]), vars)
    if k < nnodes:
        r = nodenames[k]
        for i in range(nlinks):
            vars = [n+str(i), r+'0', r+'1']
            csp.add_constraint(create_if_a_then_bc_maps_to(links[i][1]), vars)

# for i in range(nlinks):
#     lir = links[i][1]
#     for j in range(nlinks):
#         rjl = links[j][0]
#         if rjl != lir:
#             for k in range(nlinks - 1):
#                 li = linknames[k]+str(i)
#                 rj = linknames[k+1]+str(j)
#                 csp.add_constraint(not_all_2, [li,rj])

l0 = linknames[0]
r1 = linknames[nlinks-1]
csp.fix_variable(l0+str(0), 1)
csp.fix_variable(r1+str(nlinks-1), 1)
for i in range(nlinks-1):
    csp.fix_variable(l0+str(i+1), 0)
    csp.fix_variable(r1+str(i), 0)

bqm = dwavebinarycsp.stitch(csp)

print('variables:', len(csp.variables))
print('constraints:', len(csp.constraints))
print('bqm:', len(bqm))

# sampler = dimod.reference.samplers.RandomSampler()
# response = sampler.sample(bqm, num_reads=10)

sampler = dimod.reference.samplers.SimulatedAnnealingSampler()
response = sampler.sample(bqm)

# sampler = EmbeddingComposite(DWaveSampler())         # doctest: +SKIP
# response = sampler.sample(bqm, num_reads=4)         # doctest: +SKIP

sample = next(response.samples())      # doctest: +SKIP
if not csp.check(sample):              # doctest: +SKIP
    print("Failed to color map")
else:
    print("Solved!")


result = {}
for name in linknames:
    if name == l0:
        result[name] = 0
    elif name == r1:
        result[name] = nlinks - 1
    else:
        for i in range(nlinks):
            if sample[name+str(i)]:
                result[name] = i
                break
print(result)
