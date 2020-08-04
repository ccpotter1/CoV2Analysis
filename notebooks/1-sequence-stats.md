# Coronavirus Genome Calculations
    _calculates the mean and standard deviation of the lengths and gc content of your coronavirus genomes_   

I downloaded the coronavirus geneomes from NCBI
(see `data/data.md` for more details)

```julia
using BioinformaticsBISC195
genomes = parse_fasta("data/CoV_sequences.fasta") #compliles all of the genomes into 2 vectors (second vector contains sequences)
```

```julia
function seq_length(vector) #turns sequence vector into a vector containing the sequence lengths 
    length_vector = []
    for i in vector
       push!(length_vector, length(i))
    end
    return length_vector
end
```

```julia
using Statistics
meanlength= mean(seq_length(genomes[2])) #calculate mean length of the genomes
```
mean length= 29846.103770063455

```julia
using Statistics
stdlength = std(seq_length(genomes[2])) #calculate standard deviation of the genome lengths

```
standard deviation = 91.03981035873902

```julia
using BioinformaticsBISC105
function gc_contentvector(vector) #does gc_content on each part of the vector (becasue gc_content only accepts strings)
    gcvector=[]
    for i in vector
       push!(gcvector, gc_content(i)) 
    end
    return gcvector
end
```

```julia
using Statistics
mean_GC_content = mean(gc_contentvector(genomes[2])) #calculate mean gc content of all of the sequences
```
mean gc content = 0.38300315232210846

```julia
using Statistics
std_GC_content = std(gc_contentvector(genomes[2])) #calculate standard deviation of gc contents of all of the sequences
```
standard deviation of gc content = 0.009368618726296365

```julia
min_length= minimum(seq_length(genomes[2])) #find length of shortest sequence
```
minimum length = 29013

```julia
max_length = maximum(seq_length(genomes[2])) #find length of longest sequence
```
maximum length = 30484

```julia
histogram(seq_length(genomes[2])) #makes histogram of all the sequence lengths of the cov genomes
```
I first made a plot using all of the genomes and then used the remove_short_genomes function to get rid of all genomes that had less then 30000 base pairs. I then made a plot using only the genomes longer than 300000 bp.

```julia
plot!(legend=false, xaxis="Genome Length", yaxis="Number of Genomes",title="CoV Genomes") #remove the legend from the graph and creates axis labels and title
```

```julia
using BioinformaticsBISC195 
function unique_kmer_vector(vector, k) 
kmervector=[]
    for sequence in vector
        push!(kmervector, unique_kmers(sequence, k))
    end
return kmervector
end 
```