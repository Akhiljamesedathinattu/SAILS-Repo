def duplidict():
    sentence = input("Enter a sentence: ")
    words = sentence.split()
    unique_words = {}

    for word in words:
        unique_words[word] = 1
    result1 = unique_words.keys()
    print(result1)
    result = " ".join(unique_words.keys())

    print("Sentence after removing duplicates:")
    print(result)


def mymain():
    duplidict()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()