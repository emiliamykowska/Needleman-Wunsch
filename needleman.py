import seaborn as sns
import matplotlib.pyplot as plt

from needleman_utils import (
    choose_sequence_type,
    load_sequence,
    choose_algorithm_parameters,
    needleman_wunsch,
)


def main():

    seq_type = choose_sequence_type()

    first_sequence = load_sequence(
        "first", seq_type
    ).sequence  # to get a string represtnation of the sequence

    second_sequence = load_sequence("second", seq_type).sequence

    match_score, mismatch_score, gap_penalty = choose_algorithm_parameters()

    final_alignment1, final_alignment2, grid, path = needleman_wunsch(
        first_sequence, second_sequence, match_score, mismatch_score, gap_penalty
    )

    final_score = grid[-1][-1]

    matching_lengths = len(final_alignment1)  # same as len(final_alignement2)

    alignment_signs = ""

    for char1, char2 in zip(final_alignment1, final_alignment2):
        if char1 == "_" or char2 == "_":
            alignment_signs += " "
        elif char1 == char2:
            alignment_signs += "*"
        else:
            alignment_signs += "|"

    gap_percentage = alignment_signs.count(" ") / (matching_lengths) * 100
    alignment_percentage = alignment_signs.count("*") / (matching_lengths) * 100

    rows = [
        p[0] + 0.5 for p in path
    ]  # so seaborn shows lines through the centers of the cells
    cols = [p[1] + 0.5 for p in path]

    if len(first_sequence) < 100 and len(second_sequence) < 100:
        plt.figure(figsize=(16, 10))

        axes = sns.heatmap(
            grid,
            cmap="crest",
            xticklabels=["_"]
            + list(
                first_sequence
            ),  # becauase the first row and column don't correspond to chars from the sequences
            yticklabels=["_"] + list(second_sequence),
            annot=True,
        )

        axes.plot(cols, rows, color="red")
        plt.savefig("heatmap.png")
    else:
        plt.figure(figsize=(16, 10))
        plt.plot(cols, rows)
        plt.gca().invert_yaxis()  # because it was shifted
        plt.xlabel("First sequence")
        plt.ylabel("Second sequence")
        plt.savefig("heatmap.png")

    with open("needleman_report.txt", "w") as f:
        f.write(
            f"""Sequences compared:
    {first_sequence}

    {second_sequence}

    Parameters used: 
    Match score: {match_score}, 
    Mismatch score: {mismatch_score}, 
    Gap penalty: {gap_penalty}

    Alignment: 
    {final_alignment1}
    {alignment_signs}
    {final_alignment2}

    Final score: {final_score}
    Matching lengths: {matching_lengths}
    Alignment percentage: {alignment_percentage:.2f}
    Gap percentage: {gap_percentage:.2f}
    """
        )

    print("The plot and report were saved to the files")


if __name__ == "__main__":
    main()
