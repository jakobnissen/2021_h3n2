# Load in sequences from the flupipe pipeline results

##= Load in human results
dir = "results/200306"

dna_sequences = Dict{String, Vector{Tuple{Segment, LongDNASeq}}}()
aa_sequences = Dict{String, Vector{Tuple{Protein, LongAminoAcidSeq}}}()
passed_segments = Dict{String, Vector{Segment}}()

# First, we check manually all are OK. We accept all segments.
# T4794-H3N2_S22_L001 PB1 fails automated QC, but it seems alright.

for samplename in map(readdir("$dir/consensus")) do dirname
    String(first(split(dirname, '-')))
end
    passed_segments[samplename] = collect(instances(Segment))
end

for sample in readdir("$dir/consensus")
    samplename = String(first(split(sample, '-')))
    for (extension, T, SeqT, dict) in [
        ("fna", Segment, LongDNASeq, dna_sequences),
        ("faa", Protein, LongAminoAcidSeq, aa_sequences)
    ]
        open(FASTA.Reader, "$dir/consensus/$sample/consensus.$extension") do reader
            v = Tuple{T, SeqT}[]
            for record in reader
                id_fields = split(FASTA.identifier(record), '_')
                molecule_str = id_fields[end - (T === Segment ? 0 : 1)]
                molecule = molecule_str == "PB1F2" ? Proteins.PB1F2 : unwrap(parse(T, molecule_str))
                seq = FASTA.sequence(SeqT, record)
                push!(v, (molecule, seq))
            end
            dict[samplename] = v
        end
    end
end

##= Load in swine sequences
# V34866-H3N2_S23_L001 have bad PB2, PB2 and NP. All rest looks good.
dir = "results/200610"

for samplename in map(readdir("$dir/consensus")) do dirname
    String(first(split(dirname, '-')))
end
    passed = collect(instances(Segment))
    if samplename == "V34866"
        setdiff!(passed, [Segments.PB2, Segments.PB1, Segments.NP])
    end
    passed_segments[samplename] = passed
end

SAMPLES = passed_segments |> Filter(v -> !isempty(last(v))) ⨟ Map(first) |>
    x -> intersect(x, map(i -> first(split(i, '-')), readdir("$dir/consensus")))


for sample in readdir("$dir/consensus")
    samplename = String(first(split(sample, '-')))
    for (extension, T, SeqT, dict) in [
        ("fna", Segment, LongDNASeq, dna_sequences),
        ("faa", Protein, LongAminoAcidSeq, aa_sequences)
    ]
        open(FASTA.Reader, "$dir/consensus/$sample/consensus.$extension") do reader
            v = Tuple{T, SeqT}[]
            for record in reader
                id_fields = split(FASTA.identifier(record), '_')
                molecule_str = id_fields[end - (T === Segment ? 0 : 1)]
                molecule = molecule_str == "PB1F2" ? Proteins.PB1F2 : unwrap(parse(T, molecule_str))
                seq = FASTA.sequence(SeqT, record)
                segment = molecule isa Segment ? molecule : source(molecule)
                if in(segment, passed_segments[samplename])
                    push!(v, (molecule, seq))
                else
                    println("Skipping $samplename : $segment $molecule")
                end
                
            end
            dict[samplename] = sort!(v)
        end
    end
end

