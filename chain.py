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

links = [(1,2), (2,0), (0,0), (0,1)]

nlinks = len(links)
linknames = createNames('link', nlinks)
linkbits = createToupleSet(nlinks)

csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

def not_all_2(a, b):
    return not (a and b)

for name in linknames:
    variables = [name+str(i) for i in range(nlinks)]
    csp.add_constraint(linkbits, variables)

for i in range(nlinks):
    variables = [name+str(i) for name in linknames]
    csp.add_constraint(linkbits, variables)

for k in range(nlinks - 1):
    l = linknames[k]
    r = linknames[k+1]
    for i in range(nlinks):
        li = l+str(i)
        lir = links[i][1]
        for j in range(nlinks):
            rj = r+str(j)
            rjl = links[j][0]
            if not rjl == lir:
                csp.add_constraint(not_all_2, [li,rj])

bqm = dwavebinarycsp.stitch(csp)

print('variables:', len(csp.variables))
print('constraints:', len(csp.constraints))
print('bqm:', len(bqm))

l0 = linknames[0]
r1 = linknames[nlinks-1]
bqm.fix_variable(l0+str(0), 1)
bqm.fix_variable(r1+str(nlinks-1), 1)
for i in range(nlinks-1):
    bqm.fix_variable(l0+str(i+1), 0)
    bqm.fix_variable(r1+str(i), 0)

print(len(bqm))

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
    for i in range(nlinks):
        if sample[name+str(i)]:
            result[name] = i
            break
print(result)
