
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


# ============================================================
# 1. READ GENOME FILES
# ============================================================

genome1 = SeqIO.read("sequence_1.fasta", "fasta")
genome2 = SeqIO.read("sequence.fasta", "fasta")


# ============================================================
# 2. GENOME LENGTH
# ============================================================

length1 = len(genome1.seq)
length2 = len(genome2.seq)

print("Genome 1:", length1)
print("Genome 2:", length2)
print("Difference:", abs(length1 - length2))


# ============================================================
# 3. GC CONTENT
# ============================================================

def gc_content(sequence):

    sequence = sequence.upper()

    g = sequence.count("G")
    c = sequence.count("C")

    return ((g + c) / len(sequence)) * 100


gc1 = gc_content(genome1.seq)
gc2 = gc_content(genome2.seq)

print(f"GC 1: {gc1:.2f}%")
print(f"GC 2: {gc2:.2f}%")


# ============================================================
# 4. COUNT GENES
# ============================================================

def gene_count(record):

    count = 0

    for feature in record.features:

        if feature.type == "gene":
            count += 1

    return count


print("Genome 1 genes:", gene_count(genome1))
print("Genome 2 genes:", gene_count(genome2))


# ============================================================
# 5. WORKING WITH MULTIPLE FILES
# ============================================================

files = {
    "Human": "sequence.fasta",
    "Chimpanzee": "sequence_2.fasta"
}


# Empty dictionary to store sequences

sequences = {}


# Loop through the dictionary

for organism, filename in files.items():

    record = SeqIO.read(filename, "fasta")

    seq = record.seq

    sequences[organism] = seq

    print("\n" + "_" * 50)

    print("Organism:", organism)
    print("Record ID:", record.id)
    print("Description:", record.description)
    print("Length:", len(seq), "nt")
    print(f"GC content: {gc_fraction(seq) * 100:.2f}%")


# ============================================================
# 6. ACCESS INDIVIDUAL SEQUENCES FROM DICTIONARY
# ============================================================

human_seq = sequences["Human"]
chimp_seq = sequences["Chimpanzee"]


print("\n" + "_" * 50)

print("Human sequence length:", len(human_seq))
print("Chimpanzee sequence length:", len(chimp_seq))


# ============================================================
# 7. COMPARE HUMAN AND CHIMPANZEE
# ============================================================

length_difference = abs(
    len(human_seq) - len(chimp_seq)
)

print("Human-Chimpanzee length difference:",
      length_difference)


# GC comparison

human_gc = gc_fraction(human_seq) * 100
chimp_gc = gc_fraction(chimp_seq) * 100

print(f"Human GC content: {human_gc:.2f}%")
print(f"Chimpanzee GC content: {chimp_gc:.2f}%")