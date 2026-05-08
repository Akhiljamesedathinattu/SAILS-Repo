def  percentagedna():
    dna=str(input("Enter a DNA String: "))
    dna=dna.upper()
    l=len(dna)
    gcount=dna.count('G')
    ccount=dna.count('C')
    percentageG=(gcount/l)*100
    percentageC=(ccount/l)*100
    gcpercentage=((gcount+ccount)/l)*100
    print(percentageG)
    print(percentageC)
    print(gcpercentage)
def main():
    percentagedna()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')