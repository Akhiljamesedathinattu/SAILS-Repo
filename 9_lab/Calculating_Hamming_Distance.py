def hamming_distance(seq1, seq2):

    if len(seq1) != len(seq2):
        return "Sequences must be of equal length"

    count = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1

    return count



def mymain():



    # Example
    # seq1 = "GATTACA"
    # seq2 = "GACTATA"
    seq1=input("Enter first sequence: ")
    seq2=input("Enter second sequence: ")

    result = hamming_distance(seq1, seq2)

    print("Hamming Distance:", result)

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()