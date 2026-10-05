import numpy as np 

# Q2
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")
# BLOSUM62 has no entry for a gap, choose -4 as score. It is
# negative because a gap should cost points, and the alignment functions
# below ADD whatever the scoring function returns.
GAP = -4

# Score a pair of characters using BLOSUM62, with a fixed gap score.
def blosum_score(a, b):
    if a == "-" or b == "-":
        return GAP
    return blosum62[a, b]

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    n_rows = len(seq1)
    n_cols = len(seq2)

 #   dp_table[r, c]   = best score for aligning the first r letters of seq1
    #                      with the first c letters of seq2
    #   directions[r, c] = which move gave that best score (used later to
    #                      walk back and rebuild the alignment)
    dp_table = np.zeros((n_rows + 1, n_cols + 1), dtype=float)
    directions = np.zeros((n_rows + 1, n_cols + 1), dtype=int)

    # Scores for the three possible moves in the grid.
    DIAG = 1
    UP = 2
    LEFT = 3

    # First column: aligning the first r letters of seq1 against nothing.
    # Every letter has to be matched with a gap, so the score is the running
    # total of gap scores. Each of these cells was reached by moving UP.
    for r in range(1, n_rows + 1):
        dp_table[r, 0] = dp_table[r - 1, 0] + scoring_function(seq1[r - 1], "-")
        directions[r, 0] = UP

    # Initialise first row (align seq2 against gaps)
    for c in range(1, n_cols + 1):
        dp_table[0, c] = dp_table[0, c - 1] + scoring_function("-", seq2[c - 1])
        directions[0, c] = LEFT

# Fill the rest of the grid, one cell at a time, row by row.
    # The letter for row r is seq1[r - 1] and for column c is seq2[c - 1],
    # because the grid has an extra row/column at the start.   
    for r in range(1, n_rows + 1):
        for c in range(1, n_cols + 1):
            # Take the score from the cell diagonally
            # up-left and add the score for pairing these two letters.
            diag_val = dp_table[r - 1, c - 1] + scoring_function(seq1[r - 1], seq2[c - 1])

            # Take the score from the cell above and add a gap
            # penalty. seq1's letter is paired with a gap in seq2.
            up_val = dp_table[r - 1, c] + scoring_function(seq1[r - 1], "-")

            # Take the score from the cell on the left and
            # add a gap penalty. seq2's letter is paired with a gap in seq1.
            left_val = dp_table[r, c - 1] + scoring_function("-", seq2[c - 1])

            # The cell gets the best of the three options.
            best_val = max(diag_val, up_val, left_val)
            dp_table[r, c] = best_val

            # Tie-breaking order: Diag -> Up -> Left
            if best_val == diag_val:
                directions[r, c] = DIAG
            elif best_val == up_val:
                directions[r, c] = UP
            else:
                directions[r, c] = LEFT

    # Trace back from bottom right corner to origin
    aligned_a = []
    aligned_b = []
    r, c = n_rows, n_cols
    
    # Keep going until we reach the top-left corner (0, 0).
    while r > 0 or c > 0:
        step = directions[r, c]

        if step == DIAG:
            # Both letters were paired with each other. Use both, then move
            # diagonally up-left.
            aligned_a.append(seq1[r - 1])
            aligned_b.append(seq2[c - 1])
            r -= 1
            c -= 1
        elif step == UP:
            # seq1's letter was paired with a gap. Move up one row.
            aligned_a.append(seq1[r - 1])
            aligned_b.append("-")
            r -= 1
        elif step == LEFT:
            # seq2's letter was paired with a gap. Move left one column.
            aligned_a.append("-")
            aligned_b.append(seq2[c - 1])
            c -= 1
    # The lists were built from the end of the sequences to the start, so
    # reverse them and join the letters into strings.
    final_seq1 = "".join(reversed(aligned_a))
    final_seq2 = "".join(reversed(aligned_b))
    final_score = float(dp_table[n_rows, n_cols])

    return final_seq1, final_seq2, final_score
def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    # Lengths of the two sequences (rows go with seq1, columns with seq2).
    row_count = len(seq1)
    col_count = len(seq2)

    # first row and column stay 0
    dp = np.zeros((row_count + 1, col_count + 1), dtype=float)
    trace = np.zeros((row_count + 1, col_count + 1), dtype=int)

    # Scores for each move
    DIAG = 1
    UP = 2
    LEFT = 3

    # Track maximum score location for traceback start, can end in any cell 
    highest_score = 0.0
    start_r = 0
    start_c = 0

    # Fill grid
    for r in range(1, row_count + 1):
        for c in range(1, col_count + 1):
            char1 = seq1[r - 1]
            char2 = seq2[c - 1]

            match_score = dp[r - 1, c - 1] + scoring_function(char1, char2)
            gap_seq2 = dp[r - 1, c] + scoring_function(char1, "-")
            gap_seq1 = dp[r, c - 1] + scoring_function("-", char2)

            # Local alignment clips negative values to 0
            best = max(0.0, match_score, gap_seq2, gap_seq1)
            dp[r, c] = best

            #  Records which option won, with the zero check first.
            if best == 0.0:
                trace[r, c] = 0
            elif best == match_score:
                trace[r, c] = DIAG
            elif best == gap_seq2:
                trace[r, c] = UP
            else:
                trace[r, c] = LEFT

            # Update position of global maximum in the grid
            if best > highest_score:
                highest_score = best
                start_r = r
                start_c = c

    # Reconstruct alignment starting from the highest score down to 0
    align_a = []
    align_b = []
    curr_r = start_r
    curr_c = start_c

    # Stop at the edge of the grid or as soon as the score reaches 0.
    while curr_r > 0 and curr_c > 0 and dp[curr_r, curr_c] > 0:
        step = trace[curr_r, curr_c]

        if step == DIAG:
            align_a.append(seq1[curr_r - 1])
            align_b.append(seq2[curr_c - 1])
            curr_r -= 1
            curr_c -= 1
        elif step == UP:
            align_a.append(seq1[curr_r - 1])
            align_b.append("-")
            curr_r -= 1
        elif step == LEFT:
            align_a.append("-")
            align_b.append(seq2[curr_c - 1])
            curr_c -= 1
        else:
            break

    # Reverse letters and join them into strings.
    final_seq1 = "".join(reversed(align_a))
    final_seq2 = "".join(reversed(align_b))

    # The score of a local alignment is the highest value in the grid.
    return final_seq1, final_seq2, float(highest_score)

def scoring_function_simple(aa_i, aa_j):
    # Simple match/mismatch helper
    if aa_i == aa_j:
        return 1
    return -1

## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)