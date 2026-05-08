# Credit:  Introducing Python - Bill Lubanovic
# structure of a python program
import time
import math

def basics():
    # Whitespaces are used to structure the python code
    # single line comments with a # symbol

    # line continuation using \ character
    alphabet = 'abcdefg' + \
               'hijklmnop' + \
               'qrstuv' + \
               'wxyz'
    print(alphabet)
    s = 1 + 2 + \
        3
    print(s)

    ### if, elif and else
    cricket = True
    wc = True
    if cricket:
        if wc:
            print("It's cricket world cup! ")
        else:
            print("It's is not a world cup!")
    else:
        if wc:
            print("It's not a world cup of some other sport!")
        else:
            print("It's not a world cup!")

    # if more than two possibilities, use elif
    location = "bangalore"
    if location == "Chennai":
        print("The match is in Chennai")
    elif location == "Mumbai":
        print("The match is in Mumbai")
    elif location == "Delhi":
        print("The match is in Delhi")
    else:
        print("The match location is ", location)

    # comparison operators
    # equality: ==, inequality: !=, less than: <, less than or equal: <=, greater than: >
    # greater than or equal: >=, membership: in
    # All these comparison operators return boolearn value, True or False
    # use and, or, not for combining conditions
    a = 8
    if (a > 5) and (a < 10):
        print("Value of a is between 5 and 10")
    if (a > 5) and not (a == 7):
        print("Value of a is greater than 5, but not 7")

    #
    some_list = []
    if some_list:
        print("There's something in here")
    else:
        print("Hey, it's empty!")

    # do multiple comparisons with in
    letter = 'o'
    if letter == 'a' or letter == 'e' or letter == 'i' \
            or letter == 'o' or letter == 'u':
        print(letter, 'is a vowel')
    else:
        print(letter, 'is not a vowel')

    # better way
    vowels = 'aeiou'
    letter = 'o'
    if letter in vowels:
        print(letter, 'is a vowel')

    # in works with any iterables
    vowel_set = {'a', 'e', 'i', 'o', 'u'}
    print(letter in vowel_set)
    vowel_list = ['a', 'e', 'i', 'o', 'u']
    print(letter in vowel_list)
    vowel_tuple = ('a', 'e', 'i', 'o', 'u')
    print(letter in vowel_tuple)
    vowel_dict = {'a': 'apple', 'e': 'elephant', 'i': 'impala', 'o': 'ocelot', 'u': 'unicorn'}
    print(letter in vowel_dict)

    ### Repeat with while
    count = 1
    while count <= 5:
        print(count)
        count = count + 1

    # infinite loop with a break
    while True:
        stuff = input("String to capitalize [type q to quit]: ")
        if stuff == "q":
            break
        print(stuff.capitalize())

    # skip ahead with continue
    while True:
        value = input("Integer, please [q to quit]: ")
        if value == 'q':
            break
        number = int(value)
        if number % 2 == 0:  # even
            continue
        print(number, 'squared is', number * number)

    # while - else -- nonintuitive
    numbers = [1, 3, 5]
    position = 0
    while position < len(numbers):
        number = numbers[position]
        if number % 2 == 0:
            print('Found even number', number)
            break
        position += 1
    else:  # break not called
        print('No even number found')

def iterators():
    ### iterators - traversing without knowling the length
    players = ["Rohit", "Gill", "Virat", "Rahul"]
    current = 0
    while current < len(players):
        print(players[current])
        current += 1
    # here is a better pythonic way
    for x in players:
        print(x)
    # tuples, lists, strings, dict, sets are iterables, they produce an item
    # strings produce a character
    word = 'cat'
    for letter in word:
        print(letter)
    # iterating over a dictionary returns keys
    players = {'batsman': 'Rohit', 'bowler': 'bumrah', 'keeper': 'ishan'}
    for i in players:
        print(i)
    # for iterating over valueso of a dict
    for i in players.values():
        print(i)
    # return both key and value in a tuple
    for item in players.items():
        print(item)

    for type, name in players.items():
        print(type, name)

    # break and continue works with for similar to while

    # Similar to while, for has an optional else that checks if the for completed normally.
    # If break was not called, the else statement is run.

    ### iterate multiple sequences with zip
    days = ['Monday', 'Tuesday', 'Wednesday']
    fruits = ['banana', 'orange', 'peach']
    drinks = ['coffee', 'tea', 'milk']
    desserts = ['tiramisu', 'ice cream', 'pie', 'pudding']
    for day, fruit, drink, dessert in zip(days, fruits, drinks, desserts):
        print(day, ": drink", drink, "- eat", fruit, "- enjoy", dessert)
    # zip() stops when the shortest sequence is done. One of the lists (desserts) was longer than the others,
    # so no one gets any pudding unless we extend the other lists.

    # english to french dict
    english = 'Monday', 'Tuesday', 'Wednesday'
    french = 'Lundi', 'Mardi', 'Mercredi'
    print(list(zip(english, french)))
    print(dict(zip(english, french)))

    # generate number sequences with range()
    # The range() function returns a stream of numbers within a specified range.
    # You use range() similar to how to you use slices: range( start, stop, step ).
    # If you omit start, the range begins at 0. The only required value is stop;
    # as with slices, the last value created will be just before stop.
    # The default value of step is 1, but you can go backward with -1.

    # like zip, range returns an iterable object, so you need to use in to step thro' its values
    for x in range(0, 3):
        print(x)

    print(list(range(0, 3)))

    for x in range(2, -1, -1):
        print(x)

    print(list(range(2, -1, -1)))

    print(list(range(0, 11, 2)))

def comprehensions():
    # A comprehension is a compact way of creating a Python data structure from one or more iterators.
    number_list = []
    number_list.append(1)
    number_list.append(2)
    number_list.append(3)
    number_list.append(4)
    number_list.append(5)
    print(number_list)
    # another way
    number_list = []
    for number in range(1, 6):
        number_list.append(number)
    print(number_list)
    # yet another way
    number_list = list(range(1, 6))
    print(number_list)
    # Pythonic way
    # [ expression for item in iterable ]
    number_list = [number for number in range(1, 6)]
    print(number_list)
    number_list = [number - 1 for number in range(1, 6)]
    print(number_list)
    # [ expression for item in iterable if condition ]
    a_list = [number for number in range(1, 6) if number % 2 == 1]
    print(a_list)
    # nested loops
    rows = range(1, 4)
    cols = range(1, 3)
    for row in rows:
        for col in cols:
            print(row, col)
    # lets use a comprehension and assign it to the variable cells, making it a list of (row, col) tuples:
    rows = range(1, 4)
    cols = range(1, 3)
    cells = [(row, col) for row in rows for col in cols]
    for cell in cells:
        print(cell)
    for row, col in cells:
        print(row, col)

    # dictionary comprehension
    # {key_expression: value_expression for expression in iterable}
    word = 'letters'
    letter_counts = {letter: word.count(letter) for letter in word}
    print(letter_counts)
    # more pythonic - no need to repeat e and t twice
    letter_counts = {letter: word.count(letter) for letter in set(word)}

    # set comprehension
    # {expression for expression in iterable}
    a_set = {number for number in range(1, 6) if number % 3 == 1}
    print(a_set)

    # no tuple comprehension, but generator comprehensin exists
    number_thing = (number for number in range(1, 6))
    print(number_thing)
    for number in number_thing:
        print(number)
    number_list = list(number_thing)
    print(number_list)

def do_nothing():
    pass # says that the function does nothing

def make_a_sound():
    print('Quack')

def agree():
    return True

def commentary(color):
    if color == 'red':
        return "It's a tomato."
    elif color == 'green':
        return "It's a green pepper"
    elif color == 'bee purple':
        return "bee purple color"
    else:
        return "Not heard of this color " + color + '.'

def players(bat, bowl, keep, all="Hardik"):
    return {'batsman': bat, 'bowler': bowl, 'keeper': keep, 'allrounder':all}

def scope(a, result = []):
    result.append(a)
    print(result)

# When used inside the function with a parameter,
# an asterisk groups a variable number of positional arguments into a tuple of parameter values.
def print_args(*args):
    print('Positional argument tuple:', args)

# If your function has required positional arguments as well, *args goes at the end and grabs all the rest
def print_more(required1, required2, *args):
    print('Need this one:', required1)
    print('Need this one too:', required2)
    print('All the rest:', args)

# You can use two asterisks (**) to group keyword arguments into a dictionary,
# where the argument names are the keys, and their values are the corresponding dictionary values.
def print_kwargs(**kwargs):
    print('Keyword arguments:', kwargs)
    print(kwargs.keys())

def echo(anything):
    'echo returns its input argument'
    return anything

def add_args(arg1, arg2):
    print(arg1 + arg2)

def run_something_with_args(func, arg1, arg2):
    func(arg1, arg2)

def sum_args(*args):
    return sum(args)

def run_with_positional_args(func, *args):
    return func(*args)

def outer(a, b):
    def inner(c, d):
        return c + d
    return inner(a, b)

# inner2() uses the outer saying parameter directly instead of getting it as an argument.
# knights2() returns the inner2 function name instead of calling it.
def knights2(saying):
    def inner2():
        return "We are the knights who say: '%s'" % saying
    return inner2


def square_numbers(numbers, func):
    for number in numbers:
        print(func(number))

def square_function(x):
    return x*x

def functions():
    # code reuse
    # take any number of input parameters and return any number of output results
    # you can define a function and call a function
    # function is defined using def keyword, function name same as variable name rules
    # () and :
    # The next line need to be indented

    do_nothing() # calling a function
    make_a_sound()
    if agree():
        print('Splendid')
    else:
        print('That was unexpected')
    comment = commentary('blue')
    print(comment)

    # None is useful
    value = do_nothing()
    print(value)
    if value:
        print("It's some thing")
    else:
        print("It's no thing")

    # positional arguments
    # values are copied to corresponding parameters
    # IMP: you need to remember the meaning of each position
    p = players('Rohit', 'Bumrah', 'Ishan')
    print(p)

    # keyword arguments
    p = players(bat='Rohit', bowl='Bumrah', keep='Ishan')
    print(p)

    # default parameter values
    p = players(bat='Rohit', bowl='Bumrah', keep='Ishan')
    print(p)
    p = players(bat='Rohit', bowl='Bumrah', keep='Ishan', all='Jadeja')
    print(p)

    # check the scope
    scope(5)
    scope(6)  # both 5 and 6 are present in result within the function

    # gather positional arguments with *
    print_args()
    print_args(3, 2, 1, 'wait!', 'uh...')
    print_more('Rohit', 'Rahul', 'Virat', 'Gill', 'Jadeja')

    # gather keyword arguments with **
    z
    # docstrings
    help(echo)

    # Functions are first class citizens
    # Python mantra, everything is an object. This includes numbers, strings, tuples, lists,
    # dictionaries—and functions, as well. Functions are first-class citizens in Python.
    # You can assign them to variables, use them as arguments to other functions, and return them from functions.
    # This gives you the capability to do some things in Python that are
    # difficult-to-impossible to carry out in many other languages.
    type(add_args)
    run_something_with_args(add_args, 5, 9)
    run_with_positional_args(sum_args, 1, 2, 3, 4)
    # You can use functions as elements of lists, tuples, sets, and dictionaries.
    # Functions are immutable, so you can also use them as dictionary keys.

    # inner functions
    print(outer(4, 7))

    # closures
    # An inner function can act as a closure. This is a function that is dynamically generated by another function
    # and can both change and remember the values of variables that were created outside the function.
    a = knights2('Duck')
    print(type(a))
    b = knights2('Hasenpfeffer')
    print(type(b))
    # invoking functions
    print(a())
    print(b())


    # Anonymous functions - the lambda() function
    # In Python, a lambda function is an anonymous function expressed as a single statement.
    # You can use it instead of a normal tiny function.
    numbers = [1, 2, 3, 4]
    square_numbers(numbers, square_function)

    # The lambda takes one argument, which we call word here.
    # Everything between the colon and the terminating parenthesis is the definition of the function.
    square_numbers(numbers, lambda n: n * n)

    ####  properties of functions
    # A function is an instance of the Object type.
    # You can store the function in a variable
    def shout(text):
        return text.upper()
    print(shout('Hello'))
    yell = shout
    print(yell('Hello'))

    # You can pass function as a parameter to another function
    def whisper(text):
        return text.lower()

    def greet(func):
        # storing the function in a variable
        greeting = func("""Hi, I am created by a function passed as an argument.""")
        print(greeting)

    greet(shout)
    greet(whisper)

    # You can return the function from a function
    def create_adder(x):
        k = 5

        def adder(y):
            print(k)
            return x + y

        return adder

    add_15 = create_adder(15)
    print(add_15(10))

    # You can store them in data structures such as hash table, lists, etc.


def my_range(first=0, last=10, step=1):
    number = first
    while number < last:
        yield number
        number += step

def generators():
    print(sum(range(1, 101)))
    print(sum(my_range(1, 101)))
    ranger = my_range(1, 101)
    for i in ranger:
        print(i)

def decorators():

    # https://realpython.com/primer-on-python-decorators/#simple-decorators-in-python

    # ### simple decorator
    # def my_decorator(func):
    #     def my_wrapper():
    #         print("Something is happening before the function is called.")
    #         func()
    #         print("Something is happening after the function is called.")
    #
    #     return my_wrapper
    #
    # def print_hello():
    #     print("Hello!")
    #
    # # make print_hello refer to a function returned by the my_decorator function.
    # print_hello = my_decorator(print_hello)
    # # invoke print_hello() to see that the original implementation is augmented with other behaviours
    # print_hello()

    ### Pythonic pie-syntax with @ symbol for decorators

    def my_decorator(func):
        def my_wrapper():
            print("Something is happening before the function is called.")
            func()
            print("Something is happening after the function is called.")

        return my_wrapper

    @my_decorator
    def print_hello():
        print("Hello!")

    print_hello()


    # decorator to calculate duration
    # taken by any function.
    def calculate_time(func):
        # added arguments inside the inner1,
        # if function takes any arguments,
        # can be added like this.
        def inner1(*args, **kwargs):
            # storing time before function execution
            begin = time.time()

            func(*args, **kwargs)

            # storing time after function execution
            end = time.time()
            print("Total time taken in : ", func.__name__, end - begin)

        return inner1

    @calculate_time
    def fact(num):
        # sleep 2 seconds because it takes very less time
        # so that you can see the actual difference
        time.sleep(2)
        print(math.factorial(num))

    # calling the function.
    fact(10)

    def hello_decorator(func):
        def inner1(*args, **kwargs):
            print("before Execution")

            # getting the returned value
            returned_value = func(*args, **kwargs)
            print("after Execution")

            # returning the value to the original frame
            return returned_value

        return inner1

    # adding decorator to the function
    @hello_decorator
    def sum_two_numbers(a, b):
        print("Inside the function")
        return a + b

    a, b = 1, 2

    # getting the value through return of the function
    print("Sum =", sum_two_numbers(a, b))

animal = 'fruitbat'
def print_global():
    print('inside print_global:', animal)

print('at the top level:', animal)


def change_and_print_global():
    global animal
    print('inside change_and_print_global:', animal)
    animal = 'wombat'
    print('after the change:', animal)

def change_local():
    animal = 'wombat'
    print('inside change_local:', animal, id(animal))
    # Python provides two functions to access the contents of your namespaces:
    # • locals() returns a dictionary of the contents of the local namespace.
    # • globals() returns a dictionary of the contents of the global namespace.
    print('locals:', locals())
    print('globals:', globals())

def amazing():
    '''
    Documentation for this function
    Want to see it again?
    '''
    print('This function is named:', amazing.__name__)
    print('And its docstring is:', amazing.__doc__)

def namespaces():
    print_global()
    change_and_print_global()
    change_local()

    # Uses of _ and __ in Names
    # Names that begin and end with two underscores (__) are reserved for use within
    # Python, so you should not use them with your own variables.
    amazing()

def try_except():
    short_list = [1, 2, 3]
    position = 5
    try:
        print(short_list[position])
    except:
        print('Need a position between 0 and', len(short_list) - 1, ' but got',position)

    print("continue")

    # another example
    while True:
        value = input('Position [q to quit]? ')
        if value == 'q':
            break
        try:
            position = int(value)
            print(short_list[position])
        except IndexError as err:
            print('Bad index:', position)
        except Exception as other:
            print('Something else broke:', other)


    # you can define your own exceptions using class - come back to it later

def main():


    # ### basics
    # basics()
    #
    # ### Iterators
    # iterators()
    #
    # ### Comprehensions
    # comprehensions()
    #
    ### functions
    # functions()
    #

    ### generators - A generator is a Python sequence creation object. With it, you can iterate through potentially
    # huge sequences without creating and storing the entire sequence in memory at once.
    # generators()

    ### decorators - creadit: https://www.geeksforgeeks.org/decorators-in-python/
    # decorators are used to modify the behaviour of function or class. In Decorators, functions are taken as
    # the argument into another function and then called inside the wrapper function.
    # decorators()

    ### namespaces and scope
    # Python programs
    # have various namespaces—sections within which a particular name is unique and
    # unrelated to the same name in other namespaces.
    # namespaces()

    ### try and except
    # It’s good practice to add exception handling anywhere an exception might occur to let
    # the user know what is happening.
    # If an exception occurs in some function and is not caught there, it bubbles up until it is
    # caught by a matching handler in some calling function. If you don’t provide your own
    # exception handler, Python prints an error message and some information about
    # where the error occurred and then terminates the program
    try_except()



    print("End")

# Construct to not include whole program in other includes
if __name__ == "__main__":
    main()

