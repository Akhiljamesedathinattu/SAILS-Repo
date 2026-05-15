# Function to find reverse complement of a DNA sequence

def reverse_complement():
    dna = input("Enter DNA sequence: ").upper()
    print(dna)
    print(dna[::-1])
    # Dictionary for complements
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    # Find complement sequence
    combined= ""
    for b in dna:
        combined += complement[b]



# Output
    print(combined)
    print("Reverse Complement:", combined[::-1])




def mymain():
    reverse_complement()




if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()