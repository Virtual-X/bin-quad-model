import dwavebinarycsp
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx
import matplotlib.pyplot as plt
import math

import dimod

def createNames(base, count):
    return [base+str(i)+'_' for i in range(count)]

def createToupleSet(nval):
    return { ((0,)*i) + (1,) + ((0,)*(nval - i - 1)) for i in range(nval) }

links = [(1,2), (2,0), (2,0), (0,1), (1,2), (0,1)]

def log2i(n):
    return int(math.log(n, 2))

nlinks = len(links)
linknames = createNames('link', nlinks)
linkbits = log2i(nlinks - 1) + 1

nnodes = nlinks - 1
nodenames = createNames('nodes', nnodes)
nodebits = 2

csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

def not_all_2(a, b):
    return not (a and b)

def as_number(bits):
    t = 0
    f = 1
    for i in range(len(bits)):
        if bits[i]:
            t += f
        f *= 2
    return t

def get_unpadded_bits(n):
    b = []
    while n > 0:
        b.append(n & 1)
        n >>= 1
    return b

def get_bits(n, numbits=None):
    return [(n >> i) & 1 for i in range(numbits)] \
        if numbits is not None else get_unpadded_bits(n)

def create_as_number_are_smaller_than(v):
    def as_number_are_smaller_than_v(*bits):
        return as_number(bits) < v
    return as_number_are_smaller_than_v

variables = [[name+str(i) for i in range(linkbits)] for name in linknames]

if nlinks < 2 ** linkbits:
    constr = create_as_number_are_smaller_than(nlinks)
    for bits in variables:
        csp.add_constraint(constr, bits)

def second_half_is_different_from_first(*bits):
    n = len(bits) >> 1
    for i in range(n):
        if bits[i] != bits[n + i]:
            return True
    return False

for i in range(nlinks):
    for j in range(i):
        csp.add_constraint(second_half_is_different_from_first, variables[i]+variables[j])

def create_if_firsts_are_v0_then_last_two_are_v1(v0, v1):
    def if_firsts_as_number_equal_to_v0_then_bc_maps_to_v1(*bits):
        n = len(bits) - 2
        return as_number(bits[:n]) != v0 or as_number(bits[n:]) == v1
    return if_firsts_as_number_equal_to_v0_then_bc_maps_to_v1

# for k in range(nlinks):
#     n = variables[k]
#     if k > 0:
#         l = nodenames[k-1]
#         vars = n + [l+'0', l+'1']
#         for i in range(nlinks):
#             csp.add_constraint(create_if_firsts_are_v0_then_last_two_are_v1(i, links[i][0]), vars)
#     if k < nnodes:
#         r = nodenames[k]
#         vars = n + [r+'0', r+'1']
#         for i in range(nlinks):
#             csp.add_constraint(create_if_firsts_are_v0_then_last_two_are_v1(i, links[i][1]), vars)

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
bits_r2 = get_bits(nlinks-1, linkbits)
for i in range(linkbits):
    csp.fix_variable(l0+str(i), 0)
    csp.fix_variable(r1+str(i), bits_r2[i])

print('variables:', len(csp.variables))
print('constraints:', len(csp.constraints))
bqm = dwavebinarycsp.stitch(csp)
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
        result[name] = as_number([sample[name+str(i)] for i in range(linkbits)])
print(result)
