from dwave.system import DWaveSampler, EmbeddingComposite
from dimod import BinaryQuadraticModel

pumps = [0,1,2,3]
costs = [ [36,27],
         [56, 65],
         [48, 36],
         [52, 16]]
flow = [2, 7, 3, 8]
demand = 27
time = [0,1]

x = [[f'P{p}_AM', f'P{p}_PM'] for p in pumps]

bqm = BinaryQuadraticModel('BINARY')

for p in pumps:
    for t in time:
        bqm.add_variable(x[p][t], costs[p][t])

for p in pumps:
    c1 = [(x[p][t], 1) for t in time]
    bqm.add_linear_inequality_constraint(c1,
                                         lb = 1,
                                         ub = len(time),
                                         lagrange_multiplier=5,
                                         label = 'c1_pump_'+str(p))
    
for t in time:
    c2 = [(x[p][t], 1) for p in pumps]
    bqm.add_linear_inequality_constraint(c2,
                                         constant = -3,
                                         lagrange_multiplier=1,
                                         label = 'c2_time_'+str(t))
    
c3 = [(x[p][t], flow[p]) for t in time for p in pumps]
bqm.add_linear_equality_constraint(c3,
                                   constant = -demand,
                                   lagrange_multiplier=28)

sampler = EmbeddingComposite(DWaveSampler(token = 'DEV-44af657c9e8cfc328a002636f6d047d782a67e7c'))
sampleset = sampler.sample(bqm, num_reads=1000)

sample = sampleset.first.sample
total_flow = 0
total_cost= 0

print('\n\tAM\tPM')
for p in pumps:
    printout = 'P' + str(p)
    for time in range(2):
        printout += '\t' + str(sample[x[p][time]])
        total_flow += sample[x[p][time]]*flow[p]
        total_cost += sample[x[p][time]]*costs[p][time]
    print(printout)

print('Total flow:', total_flow)
print('Total cost:', total_cost)
