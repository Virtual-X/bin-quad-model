# -*- coding: utf-8 -*-

import dwavebinarycsp
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx
import matplotlib.pyplot as plt

import dimod

# Represent the map as the nodes and edges of a graph
provinces = ['AB', 'BC', 'MB', 'NB', 'NL', 'NS', 'NT', 'NU', 'ON', 'PE', 'QC', 'SK', 'YT']
neighbors = [('AB', 'BC'), ('AB', 'NT'), ('AB', 'SK'), ('BC', 'NT'), ('BC', 'YT'), ('MB', 'NU'),
             ('MB', 'ON'), ('MB', 'SK'), ('NB', 'NS'), ('NB', 'QC'), ('NL', 'QC'), ('NT', 'NU'),
             ('NT', 'SK'), ('NT', 'YT'), ('ON', 'QC')]

# Function for the constraint that two nodes with a shared edge not both select one color
def not_both_1(v, u):
    return not (v and u)

# Function that plots a returned sample
def plot_map(sample):
    G = nx.Graph()
    G.add_nodes_from(provinces)
    G.add_edges_from(neighbors)
    # Translate from binary to integer color representation
    color_map = {}
    samples = dict(sample)
    for province in provinces:
        #   for i in range(colors):
        #     if sample[province+str(i)]:
        #         color_map[province] = i
        bit0 = samples.get(province+'0', 0)
        bit1 = samples.get(province+'1', 0)
        color_map[province] = bit0 + 2 * bit1
    # Plot the sample with color-coded nodes
    node_colors = [color_map.get(node) for node in G.nodes()]
    nx.draw_circular(G, with_labels=True, node_color=node_colors, node_size=3000, cmap=plt.cm.rainbow)
    plt.show()

# at most 10011101 (157) -> (not b0) or ((not b1 and not b2) and ((not b3 or not b4 or not b5) or (not b6)
# easier (b0<<6)+(b1<<5)+...+(b5<<1)+b6 < max

# Valid configurations for the constraint that each node select a single color
one_color_configurations = {(0, 0), (0, 1), (1, 0), (1, 1)}
colors = len(one_color_configurations)

# Create a binary constraint satisfaction problem
csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

# Add constraint that each node (province) select a single color
for province in provinces:
    variables = [province+str(i) for i in range(colors)]
    #csp.add_constraint(one_color_configurations, variables)

def ab_not_equal_cd(a, b, c, d):
    return a != c or b != d

# Add constraint that each pair of nodes with a shared edge not both select one color
constrained = set()
for neighbor in neighbors:
    v, u = neighbor
    # for i in range(colors):
    #     variables = [v+str(i), u+str(i)]
    #     csp.add_constraint(not_both_1, variables)
    variables = [v+'0',v+'1', u+'0',u+'1']
    csp.add_constraint(ab_not_equal_cd, variables)

# Convert the binary constraint satisfaction problem to a binary quadratic model
bqm = dwavebinarycsp.stitch(csp)

print('variables:', len(csp.variables))
print('constraints:', len(csp.constraints))
print('bqm:', len(bqm))

# Set up a solver using the local system’s default D-Wave Cloud Client configuration file
# and sample 50 times

# bqm = dimod.BinaryQuadraticModel.

# sampler = dimod.reference.samplers.RandomSampler()
# response = sampler.sample(bqm)

sampler = dimod.reference.samplers.SimulatedAnnealingSampler()
response = sampler.sample(bqm)

# sampler = EmbeddingComposite(DWaveSampler())         # doctest: +SKIP
# response = sampler.sample(bqm, num_reads=10)         # doctest: +SKIP

sample = next(response.samples())      # doctest: +SKIP
if not csp.check(sample):              # doctest: +SKIP
    print("Failed to color map")
else:
    print("Solved!")

plot_map(sample)
