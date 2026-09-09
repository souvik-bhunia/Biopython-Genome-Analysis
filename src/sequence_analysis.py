from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd


def summarize_sequence(record):
    """Calculate sequence-level statistics for a FASTA record."""

    sequence = str(record.seq).upper()

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    n_count = sequence.count("N")

    gc_percent = (g_count + c_count) / len(sequence) * 100
    biopython_gc = gc_fraction(sequence) * 100

    valid_bases = set("ATGCN")
    unexpected_bases = set(sequence) - valid_bases

    gc_difference = abs(gc_percent - biopython_gc)

    return {
        "sequence_id": record.id,
        "length_bp": len(sequence),
        "A_count": a_count,
        "T_count": t_count,
        "G_count": g_count,
        "C_count": c_count,
        "N_count": n_count,
        "GC_percent": round(gc_percent, 2),
        "unexpected_bases": "".join(sorted(unexpected_bases)),
        "GC_verification_difference": gc_difference,
    }


# Input FASTA files
genome_files = [
    "data/Ecoli_genome.fasta",
    "data/B.subtilis.fasta",
    "data/P putida.fasta",
]


all_summaries = []

for genome_file in genome_files:

    records = SeqIO.parse(genome_file, "fasta")

    for record in records:
        summary = summarize_sequence(record)
        all_summaries.append(summary)


# Create a structured results table
summary_df = pd.DataFrame(all_summaries)


# Export results
summary_df.to_csv(
    "results/genome_summary.csv",
    index=False
)


print("Genome sequence analysis completed.")
print()
print(summary_df)
print()
print("Results saved to: results/genome_summary.csv")
