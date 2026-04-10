# import argparse

# parser = argparse.ArgumentParser(
#     prog="needleman.py", description="Needleman-Wunsch sequence alignment"
# )

# parser.add_argument("--match", type=int, required=False)
# parser.add_argument()

from needleman_utils import (
    choose_sequence_type,
    load_sequence,
    choose_algorithm_parameters,
)

seq_type = choose_sequence_type()

first_sequence = load_sequence("first", seq_type)

second_sequence = load_sequence("second", seq_type)

match_score, mismatch_score, gap_penalty = choose_algorithm_parameters()
