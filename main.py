from Bio import SeqIO
import os

SEQ_QUEUE = []
GENE_QUEUE = []

COUNT = 1
FEATURE_ID = 1

BASE = 5
MOD = 10**9 + 7
CHAR_MAP = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

MOTIFS = {
    "TATA_box": "TATAAA",
    "CAAT_box": "CCAAT",
    "GC_box": "GGGCGG",
    "Pribnow_box": "TATAAT",
    "AP1_site": "TGACTCA",
    "CRE_site": "TGACGTCA",
    "Octamer": "ATGCAAAT",
    "Kozak_sequence": "GCCACCATGG",
    "PolyA_signal": "AATAAA",
    "Shine_Dalgarno": "AGGAGG"
}

def enqueue(seq, gene):
    SEQ_QUEUE.append(seq)
    GENE_QUEUE.append(gene)

def dequeue():
    if not SEQ_QUEUE:
        return None, None
    return SEQ_QUEUE.pop(0), GENE_QUEUE.pop(0)

def gc_content(seq):
    total = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
    if total == 0:
        return 0
    return ((seq.count("G") + seq.count("C")) / total) * 100

def rabin_karp_search(seq, pattern):
    n, m = len(seq), len(pattern)
    if m > n:
        return []

    pat_hash = 0
    win_hash = 0
    high_pow = pow(BASE, m - 1, MOD)

    for i in range(m):
        if pattern[i] not in CHAR_MAP or seq[i] not in CHAR_MAP:
            return []
        pat_hash = (pat_hash * BASE + CHAR_MAP[pattern[i]]) % MOD
        win_hash = (win_hash * BASE + CHAR_MAP[seq[i]]) % MOD

    positions = []

    for i in range(n - m + 1):
        if win_hash == pat_hash and seq[i:i + m] == pattern:
            positions.append(i)

        if i < n - m:
            left = seq[i]
            right = seq[i + m]

            if left not in CHAR_MAP or right not in CHAR_MAP:
                win_hash = 0
                continue

            win_hash = (win_hash - CHAR_MAP[left] * high_pow) % MOD
            win_hash = (win_hash * BASE + CHAR_MAP[right]) % MOD

    return positions

def detect_orfs(seq):
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            if seq[i:i+3] == START_CODON:
                start = i
                i += 3
                while i < len(seq) - 2:
                    if seq[i:i+3] in STOP_CODONS:
                        end = i + 3
                        orfs.append((start, end))
                        break
                    i += 3
            i += 3
    return orfs

def write_gff(gff, gene, feature, start, end, strand, attributes):
    gff.write(
        f"{gene}\tDayhoffs_Unzipper\t{feature}\t"
        f"{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
    )

def annotations(gff, gene, seq):
    global COUNT, FEATURE_ID

    print("=" * 50)
    print(f"Sequence Annotation {COUNT}")
    print("Gene ID    :", gene)
    print("Length     :", len(seq))
    print("GC Content : {:.2f}%".format(gc_content(seq)))

    gff.write(f"##sequence-region {gene} 1 {len(seq)}\n")

    orfs = detect_orfs(seq)
    if not orfs:
        print("No ORFs found")
    else:
        print("Detected ORFs:")
        for start, end in orfs:
            print(f" Start: {start} | End: {end} | Length: {end - start}")
            write_gff(
                gff,
                gene,
                "CDS",
                start + 1,
                end,
                "+",
                f"ID=cds{FEATURE_ID}"
            )
            FEATURE_ID += 1

    print("Detected Motifs:")
    found = False
    for name, motif in MOTIFS.items():
        hits = rabin_karp_search(seq, motif)
        if hits:
            found = True
            print(f" {name} ({motif}) -> {hits}")
            for pos in hits:
                write_gff(
                    gff,
                    gene,
                    "motif",
                    pos + 1,
                    pos + len(motif),
                    "+",
                    f"ID=motif{FEATURE_ID};Name={name};Sequence={motif}"
                )
                FEATURE_ID += 1

    if not found:
        print(" No motifs found")

    print("=" * 50)
    COUNT += 1

def load_fasta(path):
    for record in SeqIO.parse(path, "fasta"):
        enqueue(str(record.seq).upper(), record.id)

def run(path):
    output = os.path.splitext(path)[0] + ".gff"
    gff = open(output, "w")
    gff.write("##gff-version 3\n")

    while True:
        seq, gene = dequeue()
        if seq is None:
            break
        annotations(gff, gene, seq)

    gff.close()

if __name__ == "__main__":
    
    print("\n Dayhoff's Unzipper")
    print("=" * 30)
    fasta_path = input("Enter FASTA file path: ")
    load_fasta(fasta_path)
    run(fasta_path)
