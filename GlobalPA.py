import sys

# Scoring scheme
MATCH = 1
MISMATCH = -1
GAP = -1

def read_sequence(file):
    """Reads a sequence from a file."""
    with open(file) as f:
        return f.readline().strip()

def align_sequences(seq1, seq2):
    """Aligns two sequences using dynamic programming."""
    m = len(seq1) + 1
    n = len(seq2) + 1

    # Initialize matrices
    dp = [[0 for _ in range(n)] for _ in range(m)]
    backtrack = [[None for _ in range(n)] for _ in range(m)]

    # Set up gap penalties
    for i in range(1, m):
        dp[i][0] = i * GAP
    for j in range(1, n):
        dp[0][j] = j * GAP

    # Fill in the rest of the matrix
    for i in range(1, m):
        for j in range(1, n):
            match_score = dp[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            gap_i = dp[i-1][j] + GAP
            gap_j = dp[i][j-1] + GAP

            dp[i][j] = max(match_score, gap_i, gap_j)
            backtrack[i][j] = "match" if match_score == dp[i][j] else ("up" if gap_i == dp[i][j] else "left")

    return dp, backtrack

def backtrack_sequences(backtrack, seq1, seq2):
    """Backtracks to find the optimal alignment."""
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if backtrack[i][j] == "match":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif backtrack[i][j] == "up":
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1

    return aligned_seq1, aligned_seq2

def main():
    """Main function to run the program."""
    # Read sequences from files
    seq1 = read_sequence(sys.argv[1])
    seq2 = read_sequence(sys.argv[2])

    # Perform alignment
    dp, backtrack = align_sequences(seq1, seq2)
    aligned_seq1, aligned_seq2 = backtrack_sequences(backtrack, seq1, seq2)

    # Print results
    print(f"Score: {dp[-1][-1]}")
    print(f"Sequence 1: {aligned_seq1}")
    print(f"Sequence 2: {aligned_seq2}")

main()
