# Here, I compare sequence similarity and identity between sequences
# and references

##= First, convert `sequence` dict to an AA dict
refasms = sequences |> Map() do (segment, subdict)
    names = map(String, collect(keys(subdict)))
    seqs = collect(values(subdict))
    refasms = InfluenzaReport.manual_check(segment, seqs)
    (segment, Dict(zip(names, refasms)))
end |> Dict{Segment, Dict{String, InfluenzaReport.ReferenceAssembly}};

# Quality Check them to be sure
reports = []
for (segment, subdict) in refasms, (name, refasm) in subdict
    push!(reports, InfluenzaReport.segment_report_lines(segment, some(refasm), Dict{String,Any}()))
end

# Check they all passed :)
for (segment, subdict) in refasms, (name, refasm) in subdict
    if is_error(refasm.passed[]) || !unwrap(refasm.passed[])
        println(segment, " ", name)
    end
end

# Create a "proteins" dict like the "sequences" one
proteins = Dict{Protein, Dict{String, LongAminoAcidSeq}}()
for (segment, subdict) in refasms, (name, refasm) in subdict
    for orfdata in refasm.orfdata
        get!(Dict{String, LongAminoAcidSeq}, proteins, orfdata.variant)[name] = orfdata.aaseq
    end
end

##= See if there are any overlap in BLAST hits
hitnames = Set(keys(sequences[Segments.PA]))
for segment in instances(Segment)
    segment == Segments.HA && continue
    intersect!(hitnames, keys(sequences[segment]))
end
println("Length of intersecting hits: ", length(hitnames))

##= Compute nucleotide and AA-level identities
# * Between our samples
# * From repr. sample to S15071
# * To A/Kansas/14/2017
# * To closest public hit

# Somewhat randomly chosen, but it has good sequence quality throughout
representative = "V34865"
@assert in(representative, SAMPLES)

identities_dna = Dict(s => Dict() for s in instances(Segment))
identities_aa = Dict(p => Dict() for p in keys(proteins))

# Add in lowest between-samples identity
for (segment, subdict) in refasms
    pass_refasms = [v for (k,v) in subdict if k in SAMPLES]
    ident_dna = 1.0
    for i in 1:length(pass_refasms)-1, j in i+1:length(pass_refasms)
        aln = pairalign(OverlapAlignment(), pass_refasms[i].asmseq, pass_refasms[j].asmseq, Influenza.DEFAULT_DNA_ALN_MODEL).aln
        ident_dna = min(ident_dna, unwrap(Influenza.alignment_identity(aln)))
    end
    identities_dna[segment]["within"] = ident_dna

    # Group orfdata by protein
    byprotein = Dict()
    for refasm in pass_refasms
        for orfdata in refasm.orfdata
            push!(get!(Vector, byprotein, orfdata.variant), orfdata)
        end
    end

    for (prot, orfdatas) in byprotein
        ident_aa = 1.0
        for i in 1:length(orfdatas)-1, j in i+1:length(orfdatas)
            aln = pairalign(OverlapAlignment(), orfdatas[i].aaseq, orfdatas[j].aaseq, Influenza.DEFAULT_AA_ALN_MODEL).aln
            ident_aa = min(ident_aa, unwrap(Influenza.alignment_identity(aln)))
        end
        identities_aa[prot]["within"] = ident_aa
    end

end

# Add dist. S15071 to repr sample
# and A/Kansas/14/2017 to refsample
for refname in ["S15071", "A/Kansas/14/2017"]
    for (segment, subdict) in refasms
        refa = subdict[representative]
        refkey = if refname == "A/Kansas/14/2017" && segment == Segments.HA
            "A/Kansas/14/2017_3C.3a"
        else
            refname
        end
        refb = subdict[refkey]
        aln = pairalign(GlobalAlignment(), refa.asmseq, refb.asmseq, Influenza.DEFAULT_DNA_ALN_MODEL).aln
        ident = unwrap(Influenza.alignment_identity(aln))
        identities_dna[segment][refname] = ident

        for (orf1, orf2) in zip(sort(refa.orfdata, by=i -> i.variant), sort(refb.orfdata, by=i -> i.variant))
            @assert orf1.variant == orf2.variant
            aln = pairalign(GlobalAlignment(), orf1.aaseq, orf2.aaseq, Influenza.DEFAULT_AA_ALN_MODEL).aln
            identities_aa[orf1.variant][refname] = unwrap(Influenza.alignment_identity(aln))
        end
    end
end

# * To Closest public one (with A/-something name)
public_dists = Dict(s => Dict() for s in instances(Segment))

for (segment, subdict) in refasms
    refa = subdict[representative]
    names = filter(collect(keys(subdict))) do name
        startswith(name, "A/")
    end
    ident = maximum(
        map(names) do name
            aln = pairalign(GlobalAlignment(), refa.asmseq, subdict[name].asmseq, Influenza.DEFAULT_DNA_ALN_MODEL).aln
            id = unwrap(Influenza.alignment_identity(aln))
            public_dists[segment][name] = id
            id
        end
    )
    identities_dna[segment]["public"] = ident

    for orfdata in refa.orfdata
        ident_aa = maximum(
            map(names) do name
                refasm = subdict[name]
                reforf = refasm.orfdata[findfirst(i -> i.variant == orfdata.variant, refasm.orfdata)]
                aln = pairalign(GlobalAlignment(), orfdata.aaseq, reforf.aaseq, Influenza.DEFAULT_AA_ALN_MODEL).aln
                unwrap(Influenza.alignment_identity(aln))
            end
        )
        identities_aa[orfdata.variant]["public"] = ident_aa
    end
end

##= Output it
columns = ["within", "S15071", "public", "A/Kansas/14/2017"]
for (dict, outfile) in [
    (identities_dna, "dna_identities.csv"),
    (identities_aa, "aa_identities.csv"),
    ]
    open("results/$outfile", "w") do file
        println(file, "Molecule,within_group,S15071,ToPublic,A/Kansas/14/2017")
        for (molecule, subdict) in sort!(collect(dict), by=first)
            print(file, molecule)
            for column in columns
                print(file, ',', round(subdict[column]*100, digits=1))
            end
            println(file)
        end
    end
end
