# Week 1 Homework 

## Q1.1 
### How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?
1Mbp x 3x = 3Mbp
3Mbp / 100bp / read = 30,000 reads

## Q1.4
### In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
~50,000 reads with 0x coverage --> 50,000/1M * 100 = 5%
About 5% of the genome has not been sequenced

### How well does this match Poisson expectations? How well does the normal distribution fit the data?
Overall, the area of the data matches decently well within the poisson expectations, slightly better compared to the fit of normal distribution. However, there are multiple regions of the histogram that stick out from the expectation data lines. This can possibly be improved by more coverage.

## Q1.5
### In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
~5 reads with 0x coverage --> 5/1M * 100 = 0.0005%
About 0.0005% of the genome has not been sequenced

### How well does this match Poisson expectations? How well does the normal distribution fit the data?
Overall, the area of the data matches better within the poisson expectations, and slightly better compared to the fit of normal distribution. Again, there are multiple regions of the histogram that stick out from the expectation data lines. This can possibly be improved by more coverage.

## Q1.6
### In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
~50 reads with 0x coverage --> 50/1M * 100 = 0.005%
About 0.005% of the genome has not been sequenced

### How well does this match Poisson expectations? How well does the normal distribution fit the data?
Overall, the area of the data matches much better compared to previous coverages, and matches the poisson/normal distributions fairly well. Although there are some portion of the histogram that stick out from the expectated distribution lines, this may be due to formatting of the bandwidth. Regardless, greater coverage provided a better fit area of distribution.



## Q2.4
### Now, use dot to produce a directed graph. Record the command you used in your READMD.md

(graphviz) cmdb@QuantBio-10 Week_1 % ./Week_1_Homework_Q2.py > ex2_digraph.dot          
#### Used to open the Python File and "export" it as a dot file
(graphviz) cmdb@QuantBio-10 Week_1 % dot -Tpng ex2_digraph.dot > ex2_digraph.png        
#### Changed the formatting of the dot file to a png file

## Q2.5
### Assume that the maximum number of occurrences of any 3-mer in the actual genome is five. Using your graph from Step 2.4, write one possible genome sequence that would produce these reads. Record your answer in your README.md.

ATTGATTCACTTATTTGATTCAT


## Q2.6
### In a few sentences, what would it take to accurately reconstruct the sequence of the genome? Record your answer in your README.md.

In order to reconstruct the sequence of the genome, there must be sufficent coverage so be able to cover the entire genome. There must also be minimal sequencing errors to propagate relatively correct genome sequences.