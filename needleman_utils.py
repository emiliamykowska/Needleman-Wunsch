from bio_classes import ProteinSequence, DNASequence
import requests


def create_sequence(sequence, seq_type, source):
    seq_type = seq_type.lower()  # for more flexibility
    if seq_type == "protein":
        return ProteinSequence(sequence, source)
    elif seq_type == "dna":
        return DNASequence(sequence, source)
    else:
        raise ValueError(
            f"Invalid sequence type '{seq_type}'. Valid types are 'protein' and 'DNA'"
        )


def fasta_to_sequence(path, seq_type, source):
    try:
        with open(path, "r") as f:
            lines = f.readlines()

            sequence = "".join(
                line.strip() for line in lines if not line.startswith(">")
            )
    except Exception as e:
        print(f"Error while reading file: {e}")
        return None

    return create_sequence(sequence, seq_type, source)


def choose_sequence_type():
    seq_type = ""

    while seq_type not in ("DNA", "protein"):
        type_input = input(
            "Which sequence type do you want to compare (enter 1 or 2)?\n"
            "1.Protein Sequence\n"
            "2.DNA sequence\n"
        ).strip()  # strip to delete trailing whitespaces

        if type_input == "1":
            seq_type = "protein"
        elif type_input == "2":
            seq_type = "DNA"
        else:
            print("Not a valid answer. Please, try again.\n")

    return seq_type


def load_sequence(seq_number, seq_type):
    while True:
        source_input = input(
            f"Choose the input source for the {seq_number} sequence (enter 1, 2 or 3):\n"
            "1. Manual sequence input\n"
            "2. Sequence from FASTA file\n"
            "3. Sequence from NCBI\n"
        ).strip()

        if source_input == "1":
            source = "Manual"
            sequence = input("Provide the sequence:\n")
            return create_sequence(sequence, seq_type, source)
        elif source_input == "2":
            source = "FASTA"
            path = input("Provide file to the path containing the FASTA format:\n")
            sequence = fasta_to_sequence(path, seq_type, source)
            if sequence is None:
                print("Please try again.")
                continue
            else:
                return sequence
        elif source_input == "3":
            source = "NCBI"
            accession_id = input("Provide accession ID:\n")
            sequence = fetch_fasta_from_ncbi(accession_id)

            if sequence is None:
                print("Please try again.")
                continue
            else:
                return create_sequence(sequence, seq_type, source)
        else:
            print("Invalid choice. Please, try again.\n")


def fetch_fasta_from_ncbi(accession_id):
    """Fetches a nucleotide sequence from the NCBI database."""
    print(f"Fetching record {accession_id} from NCBI...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession_id}&rettype=fasta&retmode=text"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        # Parsing the FASTA format (skipping the first header line)
        lines = response.text.strip().split("\n")
        sequence = "".join(lines[1:])
        return sequence
    except Exception as e:
        print(f"Error while fetching data: {e}")
        return None


def choose_algorithm_parameters():
    while True:
        choice = input(
            "Do you want to use default Needleman-Wunsch algorith parameters (match = 1, mismatch = -1, gap penalty = -2)?\n"
            "1. Yes\n"
            "2. No\n"
        )

        if choice == "1":
            return 1, -1, -2
        elif choice == "2":
            try:
                match_score = int(input("Please provide score for matching:\n").strip())

                if match_score < 0:
                    print("Match score should be >= 0\n")
                    continue

                mismatch_score = int(
                    input("Please provide score for mismtach\n").strip()
                )
                gap_penalty = int(
                    input("Please provide score for gap penalty:\n").strip()
                )

                if gap_penalty > 0:
                    print("Gap penalty should be <= 0.\n")
                    continue

                if match_score < mismatch_score or mismatch_score < gap_penalty:
                    print(
                        "Match score should not be lower than mismatch or gap penalty. Gap penalty should not be greater than mismtach score.\n"
                    )
                    continue

                return match_score, mismatch_score, gap_penalty

            except ValueError:
                print("Please enter valid integers.\n")

        else:
            print(f"Wrong answer '{choice}'. Please, try again.\n")
