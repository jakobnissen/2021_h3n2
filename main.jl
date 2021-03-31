import Conda

# Set path to be able to run MAFFT and other tools through Julia
ENV["PATH"] = Conda.BINDIR * ':' * ENV["PATH"]

using FASTX
using CodecZlib
using BioSequences
using BioAlignments
using Influenza
using InfluenzaReport
using Transducers
using ErrorTypes

maybedir(dir::AbstractString) = isdir(dir) || mkdir(dir)

# Step 1: Run pipeline  on the two used raw datasets
# This needs to be done "manually".
# You can find the pipeline and all its scripts at https://github.com/ssi-dk/flustrain#75ebe11
# The output directories should be results/200306 and results/200610, respectively.

##= Setup
# Create directories if they do not already exist
maybedir("results")
maybedir("tmp")

# Decompress every sequence into "tmp"
function mirror_decompress(target::AbstractString, src::AbstractString)
    isdir(src) || error("Source must be directory")
    run(`cp -r $src $target`)
    for (dir, _, files) in walkdir(target)
        for file in files
            endswith(file, ".gz") && run(`gunzip $("$dir/$file")`)
        end
    end
end
mirror_decompress("tmp/ref", "ref")

# Load in sequences
include("scripts/load_sequences.jl")

# Load in references
include("scripts/load_references.jl")

# Run phylogeny
maybedir("results/phylogeny")
include("scripts/phylogeny.jl")

# Calculate identities
include("scripts/identity.jl")

# Check for mutations
include("scripts/mutations.jl")
