def digest_dna(dna, site):

    # Empty list to store fragments
    fragments = []

    # Start searching from beginning
    start = 0

    # Keep checking until no site is found
    while site in dna[start:]:

        # Find position of recognition site
        pos = dna.find(site, start)

        # Take DNA before the site
        fragment = dna[start:pos]

        # Store fragment in list
        fragments.append(fragment)

        # Move start after the recognition site
        start = pos + len(site)

    # Add remaining DNA after last cut
    fragments.append(dna[start:])

    return fragments

def
    # Input DNA sequence
    dna = "AAGAATTCCGGAATTCAT"

    # Restriction enzyme recognition site
    site = "GAATTC"

    # Function call
    result = digest_dna(dna, site)

    # Print fragments
    print("DNA Fragments:")

    for part in result:
        print(part)
def mymain():
    matrix_subtraction()



if __name__ == "__main__":
    mymain()