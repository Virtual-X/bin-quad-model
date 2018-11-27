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

def createNames(base, count):
    return [base+str(i)+'_' for i in range(count)]

# Function for the constraint that two nodes with a shared edge not both select one color
def not_all_4(i, j, v, u):
    return not (i and j and v and u)

def not_all_2(i, v):
    return not (i and v)

def if_a_and_b_then_c(a, b, c):
    return not a or not b or c

def createToupleSet(nval):
    return { ((0,)*i) + (1,) + ((0,)*(nval - i - 1)) for i in range(nval) }
    result = set()
    for i in range(nval):
        value = ((0,)*i) + (1,) + ((0,)*(nval - i - 1))
        result.add(value)
    return result

# Valid configurations for the constraint that each node select a single color
nrotval = 3
rotval = createToupleSet(nrotval)
rotnames = createNames('rot', 8)

nindval = 7 # 8
indval = createToupleSet(nindval)
indnames = createNames('ind', 8)

nfaceval = 4
faceval = createToupleSet(nfaceval)
facenames = createNames('face', 12)

# Create a binary constraint satisfaction problem
csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

# Add constraint that each node (province) select a single color
for name in rotnames:
    variables = [name+str(i) for i in range(nrotval)]
    csp.add_constraint(rotval, variables)
for name in indnames:
    variables = [name+str(i) for i in range(nindval)]
    csp.add_constraint(indval, variables)
for name in facenames:
    variables = [name+str(i+1) for i in range(nfaceval)]
    csp.add_constraint(faceval, variables)

def matchingFace(face):
    return (4 - face) if (face == 1 or face == 3) else face

def faceMatch(iFace, jFace):
    return jFace == matchingFace(iFace)
#    return (jFace + iFace == 4) if (iFace == 1 or iFace == 3) else (jFace == iFace)

def face(i, iRot):
    return pieces[3 * (i % 8) + iRot % 3]

def match(i, iRot, j, jRot):
    return faceMatch(face(i, iRot), face(j, jRot))

countI = [0 for i in range(len(indnames))]

def addConstraints(i, iRot, j, jRot, constraintIndex):
    facename = facenames[constraintIndex]
    indi = indnames[i]
    v = rotnames[i]
    indj = indnames[j]
    u = rotnames[j]
    countI[i] += 1
    countI[j] += 1
    for ii in range(nindval):
        for ri in range(nrotval):
            fi = face(ii, ri + iRot)
            fj = matchingFace(face(ii, ri + jRot))
            variables = [indi+str(ii), v+str(ri), facename+str(fi)]
            csp.add_constraint(if_a_and_b_then_c, variables)
            variables = [indj+str(ii), u+str(ri), facename+str(fj)]
            csp.add_constraint(if_a_and_b_then_c, variables)

constraintIndex = 0
for j in range(8):
    if (j & 1):
        i = ((j & 2) >> 1) ^ ((j & 4) >> 2)
        addConstraints(j - 1, 2 - i, j, 1 + i, constraintIndex)
        constraintIndex += 1
    if (j & 2):
        i = (j & 1) ^ ((j & 4) >> 2)
        addConstraints(j - 2, 1 + i, j, 2 - i, constraintIndex)
        constraintIndex += 1
    if (j & 4):
        addConstraints(j - 4, 0, j, 0, constraintIndex)
        constraintIndex += 1

# Add constraint that each pair of nodes with a shared edge not both select one color
# for neighbor in neighbors:
#     v, u = neighbor
#     for i in range(colors):
#         variables = [v+str(i), u+str(i)]
#         csp.add_constraint(not_both_1, variables)

# Convert the binary constraint satisfaction problem to a binary quadratic model
print(len(csp.constraints))

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
