import dwavebinarycsp
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx
import matplotlib.pyplot as plt

import dimod

# Represent the map as the nodes and edges of a graph
faces = ['213', '314', '423', '142', '243', '124', '413', '312']
pieces = [2, 1, 3 , 3, 1, 4 , 4, 2, 3 , 1, 4, 2 , 2, 4, 3 , 2, 4, 1 , 4, 1, 3 , 3, 1, 2]
pieces = [2,1,3,  4,3,1,  4,2,3,  4,2,1,  4,3,2,  2,4,1,  1,3,4,  1,2,3]
rotnames = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']
indnames = ['i1', 'i2', 'i3', 'i4', 'i5', 'i6', 'i7']

# Function for the constraint that two nodes with a shared edge not both select one color
def not_all_4(i, j, v, u):
    return not (i and j and v and u)

def not_all_2(i, v):
    return not (i and v)

# Valid configurations for the constraint that each node select a single color
rotval = {(0, 0, 1), (0, 1, 0), (1, 0, 0)}
nrotval = len(rotval)

nindval = 7
indval = set()
for i in range(nindval):
    value = ((0,)*i) + (1,) + ((0,)*(nindval - i - 1))
    indval.add(value)

# Create a binary constraint satisfaction problem
csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

# Add constraint that each node (province) select a single color
for name in rotnames:
    variables = [name+str(i) for i in range(nrotval)]
    csp.add_constraint(rotval, variables)
for name in indnames:
    variables = [name+str(i) for i in range(nindval)]
    csp.add_constraint(indval, variables)

def faceMatch(iFace, jFace):
    return (jFace + iFace == 4) if (iFace == 1 or iFace == 3) else (jFace == iFace)

def face(i, iRot):
    return pieces[3 * (i % 8) + iRot % 3]

def match(i, iRot, j, jRot):
    return faceMatch(face(i, iRot), face(j, jRot))

def addConstraints(i, iRot, j, jRot):
    indi = indnames[i]
    v = rotnames[i]
    if j != 7:
        indj = indnames[j]
        u = rotnames[j]
    for ii in range(nindval):
        for ij in range(8):
            if (j == 7 and ij != 7) or (j != 7 and ij == 7):
                continue
            for ri in range(nrotval):
                for rj in range(nrotval):
                    if j == 7 and rj != 0:
                        continue
                    if match(ii, iRot + ri, ij, jRot + rj):
                        continue
                    if j == 7:
                        variables = [indi+str(ii), v+str(ri)]
                        csp.add_constraint(not_all_2, variables)
                    else:
                        variables = [indi+str(ii), indj+str(ij), v+str(ri), u+str(rj)]
                        csp.add_constraint(not_all_4, variables)

for j in range(8):
    if (j & 1):
        i = ((j & 2) >> 1) ^ ((j & 4) >> 2)
        addConstraints(j - 1, 2 - i, j, 1 + i)
    if (j & 2):
        i = (j & 1) ^ ((j & 4) >> 2)
        addConstraints(j - 2, 1 + i, j, 2 - i)
    if (j & 4):
        addConstraints(j - 4, 0, j, 0)

# Add constraint that each pair of nodes with a shared edge not both select one color
# for neighbor in neighbors:
#     v, u = neighbor
#     for i in range(colors):
#         variables = [v+str(i), u+str(i)]
#         csp.add_constraint(not_both_1, variables)

# Convert the binary constraint satisfaction problem to a binary quadratic model
bqm = dwavebinarycsp.stitch(csp)

# Set up a solver using the local system’s default D-Wave Cloud Client configuration file
# and sample 50 times

# bqm = dimod.BinaryQuadraticModel.

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


# Function that plots a returned sample
def plot_map(sample):
    color_map = {}
    for name in indnames:
          for i in range(nindval):
            if sample[name+str(i)]:
                color_map[name] = i
    print('ind: ' + str(color_map))
    color_map = {}
    for name in rotnames:
          for i in range(nrotval):
            if sample[name+str(i)]:
                color_map[name] = i
    print('rot: ' + str(color_map))

plot_map(sample)
