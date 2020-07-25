# Coronavirus Genome Calculations
    _calculates the mean and standard deviation of the lengths and gc content of your coronavirus genomes_   

I downloaded the coronavirus geneomes from NCBI
(see `data/data.md` for more details)

```julia
using BioinformaticsBISC195
genomes = parse_fasta("data/CoV_sequences.fasta") #compliles all of the genomes into 2 vectors (second vector contains sequences)
end
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
end
```

```julia
using Statistics
stdlength = std(seq_length(genomes[2])) #calculate standard deviation of the genome lengths
end
```

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
end
```

```julia
using Statistics
std_GC_content = std(gc_contentvector(genomes[2])) #calculate standard deviation of gc contents of all of the sequences
end
```
