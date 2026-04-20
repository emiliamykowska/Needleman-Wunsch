import pytest
from needleman_utils import needleman_wunsch


def test_needleman_wunsh():
    final_alignment1, final_alignment2, grid, _ = needleman_wunsch(
        "ACTG", "ACTT", 1, -1, -2
    )
    final_score = grid[-1][-1]
    assert final_score == 2
    assert final_alignment1 == "ACTG"
    assert final_alignment2 == "ACTT"

    final_alignment1, final_alignment2, grid, _ = needleman_wunsch(
        "AC", "ATC", 1, -1, -2
    )
    final_score = grid[-1][-1]
    assert final_score == 0
    assert final_alignment1 == "A_C"
    assert final_alignment2 == "ATC"

    final_alignment1, final_alignment2, grid, _ = needleman_wunsch(
        "LYLIFGAWAGMVGTA", "KFGAAKVV", 3, 1, 0
    )
    final_score = grid[-1][-1]
    assert final_score == 18
    assert final_alignment1 == "LYLIFGAWAGMVGTA"
    assert final_alignment2 == "___KFGA_A_KV__V"

    final_alignment1, final_alignment2, grid, _ = needleman_wunsch(
        "MKSNIQDNCQVTNPATGHLFDLNSLKNDSGYSVAYSEKGLIYIGICGGTKNCPSGVGVCFGLTKINAGSWNSQLMYVDQVLQLVYDDGAPCPSKNALKYKSVISFVCTHDSGANNKPVFVSLDKQTCTLYFSWHTPLACEKEEPRHHHHHH",
        "LSKTQQAALHH",
        1,
        -1,
        -2,
    )
    final_score = grid[-1][-1]
    assert final_score == -269


def test_needleman_wunsh_boundary():
    final_alignment1, final_alignment2, grid, _ = needleman_wunsch(
        "ACTT", "ACTT", 1, -1, -2
    )
    final_score = grid[-1][-1]
    assert final_score == 4
    assert final_alignment1 == "ACTT"
    assert final_alignment2 == "ACTT"
