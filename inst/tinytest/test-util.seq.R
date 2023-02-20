# Load default settings for CHNOSZ
reset()

info <- "count.aa() warns about unrecognized amino acids and performs substring operations"
expect_message(count.aa("ABCDEFGHIJ"), "count.aa: unrecognized letter\\(s\\) in protein sequence: B J", info = info)
myseq <- "AAAAAGGGGG"
expect_equal(count.aa(myseq, stop = 5)[, "G"], 0, info = info)
expect_equal(count.aa(myseq, start = 6)[, "A"], 0, info = info)
expect_equal(count.aa(myseq, start = 5, stop = 6)[, c("A", "G")], c(1, 1), check.attributes = FALSE, info = info)

info <- "Nucleobase sequences can be processed with count.aa(), nucleic.formula() and nucleic.complement()"
expect_message(dna <- count.aa("ABCDEFGHIJ", type = "DNA"), "count.aa: unrecognized letter\\(s\\) in DNA sequence: B D E F H I J", info = info)
expect_equal(as.numeric(dna), c(1, 1, 1, 0), info = info)
expect_equal(nucleic.formula(dna), "C14H15N13O2", info = info)
# Nucleobases can be in any order
expect_equal(nucleic.formula(dna[, 4:1, drop = FALSE]), "C14H15N13O2", info = info)
# ACG -> UGC (RNA complement)
expect_equal(nucleic.formula(nucleic.complement(dna, "RNA")), "C13H14N10O4", info = info)

info <- "count.aa() correctly processes a longer nucleobase sequence"
seq <- "ATGTCCCGTTTCTTAGTTGCATTGGTTGCCGCACTTTTAGGAGTTGCAATTGAGATGTCCCTTCTCGTTCGCGCTCAGGGGCAGCAAACCTTGCTTTTGGCTGAAGAAAGCAAGCATTTGTCGCAATTGCGTCAACTGACTTTTGAAGGCACCAATGCCGAAGCGTATTGGTCGCCTGACGGGAAATGGTTGGTCTTTCAATCCACACGCCCACCTTACAAGGCTGACCAAATCTTCATCATGAGAGCGGATGGCTCGGGAGTTCGTGTCGTCAGCACGGGCAAAGGTCGTTGCACTTGTGCCTATTTCACGCCAGATGGCAAAGGCGTTATCTTTGCTACGACCCACCTTGCTGGACCAGAACCGCCGCAAGTGCCCAAACTGGACATTCCACGCTATGTTTGGGGCGTGTTCCCAAGTTACGAACTTTACCTGCGGCGTTTGGACACGATGGAACTTATCCGCTTGACCGATAACGAAGGCTACGACGCTGAAGCGACCATTTGCTGGAAGACTGGGCGAATTGTCTTCACAAGTTACCGCAATGGCGACCTTGACCTTTACAGCATGAAATTAGACGGCAGCGATTTGAAGCGATTGACGAAAACCATCGGCTACGAGGGCGGAGCGTTCTACTCGCCCGACGGGAAGCGGATTGTCTTCCGAGCCTATTTGCCAAAGACGCCTGACGAAATTGACGAATACAAGCGGTTGCTCCAGTTAGGCGTCATAAGCCCACCAAAGATGGAGTGGGTCGTCATGGACGCCGACGGTCGCAACATGAAGCAAATC"
counts <- data.frame(A = 190, C = 203, G = 211, T = 188)
expect_equal(as.numeric(count.aa(seq, type = "DNA")), as.numeric(counts), info = info)
