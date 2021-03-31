# In this file I load reference sequences.

##= HA phylogeny
# This is special because we have way more references for HA, because
# this segment is traditionally the most sequenced

function header(rec::FASTA.Record)
    lastindex = max(last(rec.identifier), last(rec.description))
    lastindex < 2 && return nothing
    String(rec.data[2:lastindex])
end

function parse_blast_header(s::AbstractString)
    m = match(r"^EPI\d+\|[A-Z]{2}[12]?\|([^\|]+)\|EPI_ISL_\d+\|\|A_/_H3N2$", s)
    isnothing(m) && error(s)
    return m[1] # header e.g. "A/Human/Moscow/31/2020"
end

sequences = Dict(s => Dict() for s in instances(Segment))

# Both Kansas and Perth are in the picked HA, so we skip those files.
open(FASTA.Reader, "tmp/ref/picked_HA.fna") do reader
    for record in reader
        sequences[Segments.HA][FASTA.identifier(record)] = FASTA.sequence(record)
    end
end

open(FASTA.Reader, "tmp/ref/blast/HA.fna") do reader
    for record in reader
        name = parse_blast_header(header(record))
        sequences[Segments.HA][name] = FASTA.sequence(record)
    end
end

# Add our sequences
for (name, v) in dna_sequences
    seq = v[findfirst(i -> first(i) == Segments.HA, v)][2]
    sequences[Segments.HA][name] = seq
end

##= Add other segments than HA
remaining_segments = filter(!=(Segments.HA), collect(instances(Segment)))
for segment in remaining_segments

    # Decided to omit Perth - it's too far away for interest
    for filename in ["A_Kansas_14_2017.fna"]
        open(FASTA.Reader, "tmp/ref/$filename") do reader
            record = only(filter(collect(reader)) do record
                last(split(FASTA.identifier(record), '|')) == string(segment)
            end) 
            name = replace(filename[1:end-4], "_"=>"/")
            sequences[segment][name] = FASTA.sequence(record)
        end
    end

    open(FASTA.Reader, "tmp/ref/blast/$segment.fna") do reader
        for record in reader
            # header looks like: >EPI1676244|PA|A/Denmark/211/2020|EPI_ISL_407965||A_/_H3N2
            name = split(header(record), '|')[3]
            sequences[segment][name] = FASTA.sequence(record)
        end
    end

    for (name, v) in dna_sequences
        i = findfirst(i -> first(i) == segment, v)
        isnothing(i) && continue
        seq = v[i][2]
        sequences[segment][name] = seq
    end
end