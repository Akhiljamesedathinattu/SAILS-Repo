# def comprehension():
#     fruits = ['mango', 'kiwi', 'strawberry', 'guava', 'pineapple', 'mandarin orange']
#
#     vowels = "aeiou"
#
#     fruits_with_only_two_vowels = [
#         fruit for fruit in fruits
#         if sum(1 for char in fruit if char in vowels) == 2
#     ]
#
#     print(fruits_with_only_two_vowels)
# def mymain():
#     comprehension()
#
#
# if __name__ == "__main__":
#     # print("this is the beginning of my program")
#     mymain()
def comprehension():

    fruits = ['mango','kiwi','strawberry','guava','pineapple','mandarin orange']

    vowels = "aeiou"

    fruits_with_only_two_vowels = []

    for fruit in fruits:

        vowel_count = 0

        for char in fruit:

            if char in vowels:
                vowel_count += 1

        if vowel_count == 2:
            fruits_with_only_two_vowels.append(fruit)

    print(fruits_with_only_two_vowels)


def mymain():
    comprehension()


if __name__ == "__main__":
    mymain()