from Bio import SeqIO
from collections import deque
import os


class GenomeQueue:
    def __init__(self):
        self.seq_queue = deque()
        self.gene_queue = deque()

    def enqueue(self, seq, gene):
        self.seq_queue.append(seq)
        self.gene_queue.append(gene)

    def dequeue(self):
        if len(self.seq_queue) == 0:
            return None, None
        return self.seq_queue.popleft(), self.gene_queue.popleft()


class GenomeAnnotator:
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

    def __init__(self, fasta_file):
        self.queue = GenomeQueue()
        self.count = 1
        self.base = 5
        self.mod = 10**9 + 7
        self.mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
        self.start_codon = "ATG"
        self.stop_codons = {"TAA", "TAG", "TGA"}
        name = os.path.splitext(fasta_file)[0]
        self.gff = open(name + ".gff", "w")
        self.gff.write("##gff-version 3\n")
        self.fid = 1

    def load_fasta(self, path):
        for record in SeqIO.parse(path, "fasta"):
            self.queue.enqueue(str(record.seq).upper(), record.id)

    def gc_content(self, seq):
        g = seq.count("G")
        c = seq.count("C")
        a = seq.count("A")
        t = seq.count("T")
        total = a + c + g + t
        if total == 0:
            return 0
        return ((g + c) / total) * 100

    def rabin_karp_search(self, seq, motif):
        n = len(seq)
        m = len(motif)
        if m > n:
            return []

        mh = 0
        wh = 0
        h = pow(self.base, m - 1, self.mod)

        for i in range(m):
            if seq[i] not in self.mapping:
                return []
            mh = (mh * self.base + self.mapping[motif[i]]) % self.mod
            wh = (wh * self.base + self.mapping[seq[i]]) % self.mod

        pos = []

        for i in range(n - m + 1):
            if mh == wh:
                if seq[i:i + m] == motif:
                    pos.append(i)

            if i < n - m:
                left = seq[i]
                right = seq[i + m]
                if left not in self.mapping or right not in self.mapping:
                    wh = 0
                    continue
                wh = (wh - self.mapping[left] * h) * self.base
                wh = (wh + self.mapping[right]) % self.mod

        return pos

    def detect_orfs(self, seq):
        orfs = []
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                if seq[i:i+3] == self.start_codon:
                    s = i
                    i += 3
                    while i < len(seq) - 2:
                        if seq[i:i+3] in self.stop_codons:
                            e = i + 3
                            orfs.append((s, e))
                            break
                        i += 3
                i += 3
        return orfs

    def write_gff(self, gene, ftype, start, end, strand, attr):
        self.gff.write(
            gene + "\tDayhoff's Unzipper\t" + ftype + "\t" +
            str(start) + "\t" + str(end) + "\t.\t" +
            strand + "\t.\t" + attr + "\n"
        )

    def annotate(self, gene, seq):
        print("=" * 50)
        print("Sequence Annotation", self.count)
        print("Gene ID     :", gene)
        print("Length      :", len(seq))
        print("GC Content  : {:.2f}%".format(self.gc_content(seq)))
        self.gff.write(f"##sequence-region {gene} 1 {len(seq)}\n")
        orfs = self.detect_orfs(seq)
        if len(orfs) == 0:
            print("No ORFs found")
        else:
            print("Detected ORFs:")
            for s, e in orfs:
                print(" Start:", s, "| End:", e, "| Length:", e - s)
                self.write_gff(
                    gene,
                    "CDS",
                    s + 1,
                    e,
                    "+",
                    "ID=cds" + str(self.fid)
                )
                self.fid += 1

        found = False
        print("Detected Motifs:")
        for name in self.MOTIFS:
            motif = self.MOTIFS[name]
            positions = self.rabin_karp_search(seq, motif)
            if len(positions) > 0:
                found = True
                print(" ", name, "(", motif, ") ->", positions)
                for p in positions:
                    self.write_gff(
                        gene,
                        "motif",
                        p + 1,
                        p + len(motif),
                        "+",
                        "ID=motif" + str(self.fid) +
                        ";Name=" + name +
                        ";Sequence=" + motif
                    )
                    self.fid += 1

        if not found:
            print(" No motifs found")

        print("=" * 50)
        self.count += 1

    def run(self):
        while True:
            seq, gene = self.queue.dequeue()
            if seq is None:
                break
            self.annotate(gene, seq)
        self.gff.close()


if __name__ == "__main__":
    path = input("Enter FASTA file path: ")
    annotator = GenomeAnnotator(path)
    annotator.load_fasta(path)
    annotator.run()
