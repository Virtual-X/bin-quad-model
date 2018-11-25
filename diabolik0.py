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
provinces = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']
neighbors = [('AB', 'BC'), ('AB', 'NT'), ('AB', 'SK'), ('BC', 'NT'), ('BC', 'YT'), ('MB', 'NU'),
             ('MB', 'ON'), ('MB', 'SK'), ('NB', 'NS'), ('NB', 'QC'), ('NL', 'QC'), ('NT', 'NU'),
             ('NT', 'SK'), ('NT', 'YT'), ('ON', 'QC')]

# Function for the constraint that two nodes with a shared edge not both select one color
def not_both_1(v, u):
    return not (v and u)

# Valid configurations for the constraint that each node select a single color
one_color_configurations = {(0, 0, 1), (0, 1, 0), (1, 0, 0)}
colors = len(one_color_configurations)

# Create a binary constraint satisfaction problem
csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

# Add constraint that each node (province) select a single color
for province in provinces:
    variables = [province+str(i) for i in range(colors)]
    csp.add_constraint(one_color_configurations, variables)

def faceMatch(iFace, jFace):
    return (jFace + iFace == 4) if (iFace == 1 or iFace == 3) else (jFace == iFace)

def face(i, iRot):
    return pieces[3 * i + iRot % 3]

def match(i, iRot, j, jRot):
    return faceMatch(face(i, iRot), face(j, jRot))

def addConstraints(i, iRot, j, jRot):
    v = provinces[i]
    u = provinces[j]
    for ri in range(colors):
        for rj in range(colors):
            if match(i, iRot + ri, j, jRot + rj):
                continue
            variables = [v+str(ri), u+str(rj)]
            csp.add_constraint(not_both_1, variables)

for j in range(len(provinces)):
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

# sampler = dimod.reference.samplers.SimulatedAnnealingSampler()
# response = sampler.sample(bqm)

sampler = EmbeddingComposite(DWaveSampler())         # doctest: +SKIP
response = sampler.sample(bqm, num_reads=4)         # doctest: +SKIP

sample = next(response.samples())      # doctest: +SKIP
if not csp.check(sample):              # doctest: +SKIP
    print("Failed to color map")
else:
    print("Solved!")


# Function that plots a returned sample
def plot_map(sample):
    color_map = {}
    for province in provinces:
          for i in range(colors):
            if sample[province+str(i)]:
                color_map[province] = i
    print(color_map)
    return
    G = nx.Graph()
    G.add_nodes_from(provinces)
    G.add_edges_from(neighbors)
    # Translate from binary to integer color representation
    color_map = {}
    for province in provinces:
          for i in range(colors):
            if sample[province+str(i)]:
                color_map[province] = i
    # Plot the sample with color-coded nodes
    node_colors = [color_map.get(node) for node in G.nodes()]
    nx.draw_circular(G, with_labels=True, node_color=node_colors, node_size=3000, cmap=plt.cm.rainbow)
    plt.show()

plot_map(sample)
