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
            path = input("Provide path to the file containing the FASTA format:\n")
            sequence = fasta_to_sequence(path, seq_type, source)
            if sequence is None:
                print("Please try again")
                continue
            else:
                return sequence
        elif source_input == "3":
            source = "NCBI"
            accession_id = input("Provide accession ID:\n")
            sequence = fetch_fasta_from_ncbi(accession_id, seq_type)

            if sequence is None:
                print("Please try again")
                continue
            else:
                return create_sequence(sequence, seq_type, source)
        else:
            print("Invalid choice. Please, try again\n")


def fetch_fasta_from_ncbi(accession_id, seq_type):
    """Fetches a nucleotide sequence from the NCBI database."""
    print(f"Fetching record {accession_id} from NCBI...")

    database = (
        "protein" if seq_type.lower() == "protein" else "nucleotide"
    )  # to enable to provide id to both proteins and nucleotides
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database}&id={accession_id}&rettype=fasta&retmode=text"

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
            "Do you want to use default Needleman-Wunsch algorithm parameters (match = 1, mismatch = -1, gap penalty = -2)?\n"
            "1. Yes\n"
            "2. No\n"
        )

        if choice == "1":
            return 1, -1, -2
        elif choice == "2":
            try:
                match_score = int(input("Please provide score for matching:\n").strip())

                if match_score < 0:
                    print("Match score should be >= 0\nPlease, try again.\n")
                    continue

                mismatch_score = int(
                    input("Please provide score for mismtach\n").strip()
                )
                gap_penalty = int(
                    input("Please provide score for gap penalty:\n").strip()
                )

                if gap_penalty > 0:
                    print("Gap penalty should be <= 0.\nPlease, try again.")
                    continue

                if match_score < mismatch_score or mismatch_score < gap_penalty:
                    print(
                        "Match score should not be lower than mismatch or gap penalty. Gap penalty should not be greater than mismtach score.\nPlease, try again.\n"
                    )
                    continue

                return match_score, mismatch_score, gap_penalty

            except ValueError:
                print("Please enter valid integers.\n")

        else:
            print(f"Wrong answer '{choice}'. Please, try again.\n")


def needleman_wunsch(sequence1, sequence2, match_score, mismatch_score, gap_penalty):
    # Add + 1 to lengths because the first row and column are empty
    len1 = len(sequence1) + 1
    len2 = len(sequence2) + 1

    # Create a grid filled with zeros
    grid = [[0 for _ in range(len1)] for _ in range(len2)]

    # Fill the first column and row with gap penalties
    for row_index in range(1, len2):
        grid[row_index][0] = row_index * gap_penalty
    for col_index in range(1, len1):
        grid[0][col_index] = col_index * gap_penalty

    for row_index in range(1, len2):
        for col_index in range(1, len1):
            # Going on the diagonal
            if sequence2[row_index - 1] == sequence1[col_index - 1]:
                diagonal = grid[row_index - 1][col_index - 1] + match_score
            else:
                diagonal = grid[row_index - 1][col_index - 1] + mismatch_score

            # Going down
            down = grid[row_index - 1][col_index] + gap_penalty

            # Going right
            right = grid[row_index][col_index - 1] + gap_penalty

            max_score = max(diagonal, down, right)
            grid[row_index][col_index] = max_score

    alignment1 = []
    alignment2 = []
    path = []

    # Right and bottom corner
    current_row = len2 - 1
    current_col = len1 - 1

    while current_row > 0 or current_col > 0:
        # Append each cell to the path
        path.append((current_row, current_col))

        # Recalculating the results to see which path was used to get here
        # Check if the best score was diagonal
        if current_row > 0 and current_col > 0:
            if sequence2[current_row - 1] == sequence1[current_col - 1]:
                score_diff = match_score
            else:
                score_diff = mismatch_score

            if (
                grid[current_row][current_col]
                == grid[current_row - 1][current_col - 1] + score_diff
            ):
                alignment1.append(sequence1[current_col - 1])
                alignment2.append(sequence2[current_row - 1])
                current_row = current_row - 1
                current_col = current_col - 1
                continue  # Skip the other checks and start the next loop iteration

        # Checking if the best score was above
        if (
            current_row > 0
            and grid[current_row][current_col]
            == grid[current_row - 1][current_col] + gap_penalty
        ):
            alignment1.append("_")
            alignment2.append(sequence2[current_row - 1])
            current_row = current_row - 1
        else:  # The best score was to the left
            alignment1.append(sequence1[current_col - 1])
            alignment2.append("_")
            current_col = current_col - 1

    # Flipping the lists as they are backwards after coming back
    final_alignment1 = "".join(alignment1[::-1])
    final_alignment2 = "".join(alignment2[::-1])

    path.append((0, 0))

    return (
        final_alignment1,
        final_alignment2,
        grid,
        path,
    )
