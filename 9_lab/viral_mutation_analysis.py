def mutation():
    # Viral mutations in two strains
    strain_a = {"C241T", "C3037T", "A23403G"}
    strain_b = {"C241T", "G25563T", "C3037T", "T28144C"}

    # 1. Mutations present in both strains (Intersection)
    shared_mutations = strain_a & strain_b

    # 2. Mutations only in strain_b (Difference)
    unique_to_b = strain_b - strain_a

    # 3. All unique mutations from both strains (Union)
    all_mutations = strain_a | strain_b

    # Display results
    print("Shared mutations:", shared_mutations)
    print("Unique to strain_b:", unique_to_b)
    print("All unique mutations:", all_mutations)























def mymain():
    mutation()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()