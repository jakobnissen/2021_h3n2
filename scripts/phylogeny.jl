# In this file, I make phylogenetic trees for all 8 segments.


# LOAD IN SEQUENCES

# CAT TOGETHER
open(FASTA.Writer, "tmp/HA.cat.fna") do writer
    foreach(sequences[Segments.HA]) do (name, seq)
        write(writer, FASTA.Record(name, seq))
    end
end

# Align
run(pipeline(`linsi --thread 8 tmp/HA.cat.fna`, stdout="tmp/HA.aln.fna"))

# Trim
run(`trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in tmp/HA.aln.fna -fasta -out tmp/HA.trim.fna`)

# IQ-TREE
# Revisit the model later.
run(`iqtree -s tmp/HA.trim.fna -pre results/phylogeny/HA.iqtree -T 2`)

for segment in remaining_segments
    open(FASTA.Writer, "tmp/$segment.cat.fna") do writer
        foreach(sequences[segment]) do (name, seq)
            write(writer, FASTA.Record(name, seq))
        end
    end
end

# Align
for segment in remaining_segments
    run(pipeline(`linsi --thread 8 tmp/$segment.cat.fna`, stdout="tmp/$segment.aln.fna"))
end

# Trim
for segment in remaining_segments
    run(`trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in tmp/$segment.aln.fna -fasta -out tmp/$segment.trim.fna`)
end

# IQ-TREE
# Revisit the model later.
for segment in remaining_segments
    run(`iqtree -s tmp/$segment.trim.fna -pre results/phylogeny/$segment.iqtree -T 2`)
end
