#!/usr/bin/env python3
import sys
import numpy
import scipy

# For 3x coverage
# Parameters:
genome_size = 1000000    # 1Mbp genome size
read_size = 100        #100bp reads
coverage = 3            # 3x coverage

# Pseudo Code
num_reads = int((genome_size * coverage)/read_size)

#num_reads = 30000

## use an array to keep track of the coverage at each position in the genome

genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads):
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1

numpy.savetxt("genome_coverage_3x.txt", genome_coverage)

# ## get the range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage + 1))

# ## Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# ## Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))


# For 10x coverage
# Parameters:
genome_size = 1000000    # 1Mbp genome size
read_size = 100        #100bp reads
coverage = 10            # 10x coverage

# Pseudo Code
num_reads = int((genome_size * coverage)/read_size)

#num_reads = 100000

## use an array to keep track of the coverage at each position in the genome

genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads):
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1

#numpy.savetxt("genome_coverage_10x.txt", genome_coverage)

# ## get the range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage + 1))

# ## Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# ## Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))



# For 30x coverage
# Parameters:
genome_size = 1000000    # 1Mbp genome size
read_size = 100        #100bp reads
coverage = 30            # 30x coverage

# Pseudo Code
num_reads = int((genome_size * coverage)/read_size)

#num_reads = 300000

## use an array to keep track of the coverage at each position in the genome

genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads):
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1

numpy.savetxt("genome_coverage_30x.txt", genome_coverage)

# ## get the range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage + 1))

# ## Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# ## Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))