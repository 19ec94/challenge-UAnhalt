from sympy import symbols, expand, Poly
from fractions import Fraction

def create_cycle_index_polynomial(c):
    # Create symbols for colours: c1, c2, ..., c_c
    colours = symbols(f'c1:{c+1}')  # c1 to c_c

    # Build a_i = sum of c_j^i for j=1..c, i=1..4
    a = [0]  # dummy zero index needed for index and cycle length matching
    for i in range(1, 5):
        a_i = sum(colour**i for colour in colours)
        a.append(a_i)

    # Cube rotation group cycle index polynomial
    Z = Fraction(1,24) * (
        a[1]**6 +          # identity
        3 * a[1]**2 * a[2]**2 +  # 180° edge rotations
        6 * a[1]**2 * a[4] +     # 120° corner rotations
        6 * a[2]**3 +            # 180° face rotations
        8 * a[3]**2              # 120° edge rotations
    )

    # Fully expand polynomial
    Z_expanded = expand(Z)

    # Convert to polynomial in colour variables
    Z_in_colour_variables = Poly(Z_expanded, *colours)
    return Z_in_colour_variables 

def count_cube_colourings(Z_in_colour_variables, m):
    total_distinct_colourings = 0
    relevant_colour_combinations = []
    # poly_Z.terms() returns list of (exponents_tuple, coefficient)
    for exponents, coeff in Z_in_colour_variables.terms():
        # exponents is a tuple like (e1, e2, ..., e_c)
        if sum(exponents) == 6 and all(e <= m for e in exponents):
            relevant_colour_combinations.append([exponents, coeff])
            total_distinct_colourings += coeff

    return [relevant_colour_combinations, total_distinct_colourings]

def get_integer_input(prompt, min_value=None, max_value=None):
    """
    Prompt user for an integer input within optional min and max bounds.
    Repeats until valid input is given.

    :param prompt: The input prompt string.
    :param min_value: Minimum acceptable integer (inclusive).
    :param max_value: Maximum acceptable integer (inclusive).
    :return: Validated integer input.
    """
    while True:
        try:
            value = int(input(prompt))
            if min_value is not None and value < min_value:
                print(f"Value must be at least {min_value}. Please try again.\n")
                continue
            if max_value is not None and value > max_value:
                print(f"Value must be at most {max_value}. Please try again.\n")
                continue
            return value
        except ValueError:
            print("Invalid input! Please enter a valid integer.\n")


# Usage:
print("Enter the number of colours (at least 1) to paint the cube's facets:")
c = get_integer_input("Number of colours (c) = ", min_value=1)

print("\nEnter the maximum number of repetitions allowed per colour (1 to 6):")
m = get_integer_input("Maximum repetitions per colour (m) = ", min_value=1, max_value=6)

Z_in_colour_variables = create_cycle_index_polynomial(c)
relevant_colour_combinations, total_distinct_colourings = count_cube_colourings(Z_in_colour_variables,m)

for colour_combination, coefficient in relevant_colour_combinations:
    combination_str = ", ".join(f"c{i+1}:{value}" for i, value in enumerate(colour_combination))
    print(f"Using the colour combination ({combination_str}), there are {coefficient} distinct ways to colour the cube.")

print(f"Total number of distinct colourings with c={c} colours and max repetitions per colour m={m} is: {total_distinct_colourings}")

