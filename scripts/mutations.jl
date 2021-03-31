# Mutations in V34865 from Flusurver
# M2 31N adamantine resistance
# NA 155H oseltamivir resistance
# NA 149A zanamivir resistance
# Pos 203+206 are involved in a-2,6/a-2,3 binding, but not the same aminoacids

##= Mutations relative to human strain S15071
# Are they all same length (yes!). That makes things easier
for (prot, subdict) in proteins
    seqa = subdict["S15071"]
    @assert all(length(seqa) == length(subdict[s]) for s in keys(subdict) if startswith(s, "V"))
end

for (prot, subdict) in sort(collect(proteins), by=first)
    seqa = subdict["S15071"]
    seqb = subdict["A/Kansas/14/2017" * (prot == Proteins.HA ? "_3C.3a" : "")]
    M = Matrix{AminoAcid}(undef, length(SAMPLES) + 2, length(seqa))
    M[1, :] .= seqa
    M[2, :] .= seqb[1:length(seqa)]
    for (i, sample) in enumerate(SAMPLES)
        if haskey(subdict, sample)
            M[i+2, :] .= subdict[sample]
        else
            M[i+2, :] .= AA_Gap
        end
    end
    for i in 1:size(M, 2)
        if length(delete!(Set(M[:,i]), AA_Gap)) != 1
            println(rpad(prot, 3), " ", rpad(i, 3), " ", join(M[:,i], ' '))
        end
    end
end

#= From above, we get these conserved mutations relative to S15071
PB2: I316M (also M in A/Kansas), R508C
HA I140S (Also S in A/Kansas), Q517H
NS1: I50L (Also L in A/Kansas)
=#