class BiologicalSequence:
    def __init__(self, sequence: str, seq_type: str, source: str):
        self.sequence = sequence.upper()  # so provided sequence is case insensitive
        self.type = seq_type
        self.source = source


class ProteinSequence(BiologicalSequence):
    def __init__(
        self, sequence: str, source: str, valid_amino_acids="ACDEFGHIKLMNPQRSTVWY"
    ):
        super().__init__(sequence, "protein", source)
        self._valid_amino_acids = valid_amino_acids

        for amino_acid in self.sequence:
            if amino_acid not in self._valid_amino_acids:
                raise ValueError(
                    f"Invalid amino acid '{amino_acid}' found in protein sequence."
                )


class DNASequence(BiologicalSequence):
    def __init__(self, sequence: str, source: str, valid_nucleotides="ACGT"):
        super().__init__(sequence, "DNA", source)
        self._valid_nucleotides = valid_nucleotides

        for nucleotide in self.sequence:
            if nucleotide not in self._valid_nucleotides:
                raise ValueError(
                    f"Invalid nucleotide '{nucleotide}' found in DNA sequence."
                )
