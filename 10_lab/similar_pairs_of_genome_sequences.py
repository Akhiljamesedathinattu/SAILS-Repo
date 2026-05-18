def comprehension():

    org1 = ["ACGTTTCA", "AGGCCTTA", "AAAACCTG"]

    org2 = ["AGCTTTGA", "GCCGGAAT", "GCTACTGA"]

    threshold = 4

    similar_pairs = []

    for seq1 in org1:

        for seq2 in org2:

            count = 0

            for a, b in zip(seq1, seq2):

                if a == b:
                    count += 1

            if count > threshold:

                similar_pairs.append((seq1, seq2))

    print(similar_pairs)


def mymain():

    comprehension()


if __name__ == "__main__":

    mymain()