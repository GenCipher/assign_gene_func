import numpy as np 

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
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)