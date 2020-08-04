# Coronavirus Genome Calculations
    _calculates the mean and standard deviation of the lengths and gc content of your coronavirus genomes_   

I downloaded the coronavirus geneomes from NCBI
(see `data/data.md` for more details)

First I made a vector with two indices: genome headers (genomes[1]) and genome sequences (genomes[2])
```julia
using BioinformaticsBISC195
genomes = parse_fasta("data/CoV_sequences.fasta") #compliles all of the genomes into 2 vectors (second vector contains sequences)
```
I then created a fucntion that finds the length of each genome

```julia
function seq_length(vector) #turns sequence vector into a vector containing the sequence lengths 
    length_vector = []
    for i in vector
       push!(length_vector, length(i))
    end
    return length_vector
end
```
I did statistics on the genome lengths 

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
found out more about the genomes by doing statistics on the GC content  
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

after calculating the minimum and maximum genome lengths, i plotted all of the genome lengths using a histogram plot

```julia
using Plots
histogram(seq_length(genomes[2])) #makes histogram of all the sequence lengths of the cov genomes
```

```julia
plot!(legend=false, xaxis="Genome Length", yaxis="Number of Genomes",title="CoV Genomes") #remove the legend from the graph and creates axis labels and title
```
I first made a plot using all of the genomes and then used the remove_short_genomes function to get rid of all genomes that had less then 29500 base pairs. I then made a plot using only the genomes longer than 300000 bp.

Since all of the genomes were longer than 29,000 bp, I removed all genomes that were shorter than 29,500 bp. 

```julia
function myisless(x)
    return length(x)< 29500 
end
```
```julia
shortgenomes = findall(myisless, genomes[2]) #creates array with the indicies of genomes that are less than 30,000 bp
deleteat!(genomes[2], shortgenomes) #removes the short genomes from the seq vector
deleteat!(genomes[1], shortgenomes) #removes the short genomes from the header vector 
```
now the sequences in genomes[2] all have at least 29,500 bp

I ran <histogram(seq_length(genomes[2]))>  and <plot!(legend=false, xaxis="Genome Length", yaxis="Number of Genomes",title="CoV Genomes")> again to get a new plot that didn't include the short genomes

I made the function unique_kmer_vector so that the unique_kmers function in Bioinformatics BISC195 can work on vectors since genomes[2] is a vector

```julia
using BioinformaticsBISC195 
function unique_kmer_vector(vector, k) 
kmervector=[]
    for sequence in vector
        append!(kmervector, unique_kmers(sequence, k)) #unique kmers returns a vector called kmers of all the kmers so append! takes the contents of the kmers vector and puts them in the kmervecotr vector
    end
return unique(kmervector)
end 
```
I ran ```Julia unique_kmer_vector(genomes[2], 3) ``` to get a vector with all of the unique kmers of length 3 

To test my kmerset_distance function that I added to BioinformaticsBISC195, I created two Sets of kmers

```julia
kmerset1 = Set(["ATTG", "GTCC", "ATCG", "CGGT"])
kmerset2 = Set(["ATTT", "GTAC", "ATAG", "CGGT"])
```

the distance between these two kmer sets is 0.8571428571428572