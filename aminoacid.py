def aminoacid():
        s1 = input("enter the DNA sequence:")
        s1 = s1.upper()
        n = len(s1)
        s2 = s1.replace('T', 'U')
        v1=n//3
        print(s2)
        print("total number of aminoacid for this sequence will be",v1AU,"and the codons and aminoacids are")
        i=0
        while(i<=n-3):
            v = s2[i:i+3]
            print(v)
            i = i + 3
            if v in['UUU','UUC']:
                print("Phenylalanine (Phe, F)")
            elif v in['UUA' ,'UUG' ,'CUU', 'CUC','CUA','CUG']:
                print("Leucine (Leu)")
            elif v in['AUU', 'AUC', 'AUA']:
                print("Isoleucine (Ile)")
            elif v in['AUG']:
                print('Methionine (Met) — Start Codon')
            elif v in['AUU', 'AUC', 'AUA']:
                print("Isoleucine (Ile)")
            elif v in['GUU', 'GUC', 'GUA', 'GUG']:
                print("Valine (Val)")
            elif v in['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']:
                print("Serine (Ser)")
            elif v in['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']:
                print("Arginine (Arg)")
            elif v in['CCU', 'CCC', 'CCA', 'CCG']:
                print("Proline (Pro)")
            elif v in['ACU', 'ACC', 'ACA', 'ACG']:
                print("Threonine (Thr)")
            elif v in['GCU', 'GCC', 'GCA', 'GCG']:
                print("Alanine (Ala)")
            elif v in['UAU', 'UAC']:
                print("Tyrosine (Tyr)")
            elif v in['CAU','CAC']:
                print("Histidine (His)")
            elif v in['CAA', 'CAG']:
                print("Glutamine (Gln)")
            elif v in['AAU', 'AAC']:
                print("Asparagine (Asn)")
            elif v in['AAA', 'AAG']:
                print("Lysine (Lys)")
            elif v in['GAU', 'GAC']:
                print("Aspartic Acid (Asp)")
            elif v in['GAA', 'GAG']:
                print("Glutamic Acid (Glu)")
            elif v in['UGU', 'UGC']:
                print("Cysteine (Cys)")
            elif v in['UGG']:
                print("Tryptophan (Trp)")
            elif v in['GGU', 'GGC', 'GGA', 'GGG']:
                print("Glycine (Gly)")
            elif v in['UAA','UAG','UGA']:
                print("Stop Codons")


            
   
   
   
   
def mymain():
    aminoacid()


if __name__ == "__main__":
    mymain()



