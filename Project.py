from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# ============================================================
# 1. FILES AND ORGANISMS
# ============================================================

files = {
    "Genome 1": "sequence_1.fasta",
    "Genome 2": "sequence.fasta",
    "Genome 3": "sequence_2.fasta"
}


# ============================================================
# 2. READ FASTA FILES
# ============================================================

genomes = {}

for organism, filename in files.items():

    record = SeqIO.read(filename, "fasta")

    genomes[organism] = record


# ============================================================
# 3. DISPLAY BASIC SEQUENCE INFORMATION
# ============================================================

print("=" * 60)
print("              GENOME ANALYSIS")
print("=" * 60)

for organism, record in genomes.items():

    sequence = record.seq

    length = len(sequence)
    gc = gc_fraction(sequence) * 100

    print("\n" + "-" * 60)
    print("Organism:", organism)
    print("File:", files[organism])
    print("Record ID:", record.id)
    print("Description:", record.description)
    print("Length:", length, "nt")
    print(f"GC Content: {gc:.2f}%")


# ============================================================
# 4. NUCLEOTIDE COMPOSITION
# ============================================================

def nucleotide_composition(sequence):

    sequence = str(sequence).upper()

    a = sequence.count("A")
    t = sequence.count("T")
    g = sequence.count("G")
    c = sequence.count("C")

    return a, t, g, c


print("\n" + "=" * 60)
print("              NUCLEOTIDE COMPOSITION")
print("=" * 60)

for organism, record in genomes.items():

    a, t, g, c = nucleotide_composition(record.seq)

    print("\n" + organism)
    print("A:", a)
    print("T:", t)
    print("G:", g)
    print("C:", c)


# ============================================================
# 5. COMPARE GENOME LENGTHS
# ============================================================

lengths = {}

for organism, record in genomes.items():

    lengths[organism] = len(record.seq)


print("\n" + "=" * 60)
print("              GENOME LENGTH COMPARISON")
print("=" * 60)

for organism, length in lengths.items():

    print(organism, ":", length, "nt")


# Difference between Genome 1 and Genome 2

difference_1_2 = abs(
    lengths["Genome 1"] - lengths["Genome 2"]
)

print("\nDifference between Genome 1 and Genome 2:",
      difference_1_2, "nt")


# Difference between Genome 1 and Genome 3

difference_1_3 = abs(
    lengths["Genome 1"] - lengths["Genome 3"]
)

print("Difference between Genome 1 and Genome 3:",
      difference_1_3, "nt")


# Difference between Genome 2 and Genome 3

difference_2_3 = abs(
    lengths["Genome 2"] - lengths["Genome 3"]
)

print("Difference between Genome 2 and Genome 3:",
      difference_2_3, "nt")


# ============================================================
# 6. GC CONTENT COMPARISON
# ============================================================

gc_values = {}

for organism, record in genomes.items():

    gc_values[organism] = gc_fraction(record.seq) * 100


print("\n" + "=" * 60)
print("              GC CONTENT COMPARISON")
print("=" * 60)

for organism, gc in gc_values.items():

    print(f"{organism}: {gc:.2f}%")


# ============================================================
# 7. BASIC SEQUENCE VALIDATION
# ============================================================

def check_sequence(sequence):

    sequence = str(sequence).upper()

    valid_bases = set("ATGC")

    invalid_bases = set(sequence) - valid_bases

    return invalid_bases


print("\n" + "=" * 60)
print("              SEQUENCE VALIDATION")
print("=" * 60)

for organism, record in genomes.items():

    invalid = check_sequence(record.seq)

    if len(invalid) == 0:

        print(organism, ": Valid DNA sequence")

    else:

        print(
            organism,
            ": Invalid/ambiguous bases found:",
            invalid
        )


# ============================================================
# 8. SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("                    SUMMARY")
print("=" * 60)

for organism, record in genomes.items():

    sequence = record.seq

    print(
        f"{organism}: "
        f"{len(sequence)} nt | "
        f"GC: {gc_fraction(sequence) * 100:.2f}%"
    )

print("\nAnalysis completed successfully.")