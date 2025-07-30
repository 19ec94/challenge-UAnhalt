from sympy import symbols, expand, Poly
from fractions import Fraction

def expand_polynomials(c):
    # Create symbols for colors: c1, c2, ..., c_c
    colors = symbols(f'c1:{c+1}')  # c1 to c_c

    # Build a_i = sum of c_j^i for j=1..c, i=1..4
    a = [0]  # dummy zero index
    for i in range(1, 5):
        a_i = sum(color**i for color in colors)
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

    # Convert to polynomial in color variables
    poly_Z = Poly(Z_expanded, *colors)
    return poly_Z

def count_cube_colorings(poly_Z, m):
    total = 0
    # poly_Z.terms() returns list of (exponents_tuple, coefficient)
    for exponents, coeff in poly_Z.terms():
        # exponents is a tuple like (e1, e2, ..., e_c)
        if sum(exponents) == 6 and all(e <= m for e in exponents):
            print(exponents, coeff)
            total += coeff

    return total

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


# Usage with clear prompts and validation:
print("Enter the number of colours (at least 1) to paint the cube's facets:")
c = get_integer_input("Number of colours (c) = ", min_value=1)

print("\nEnter the maximum number of repetitions allowed per colour (1 to 6):")
m = get_integer_input("Maximum repetitions per colour (m) = ", min_value=1, max_value=6)

# Input:
#while True:
#    try: 
#        print(f"Enter the number of colours (1 <= c) to paint the cube's facets: ")
#        c = int(input("c = "))
#        if c < 1:
#            print("The number of colours must be at least greater than 1.\n")
#            continue
#
#        break
#
#    except ValueError:
#        print("Please enter a valid integer")
#        
#
#while True:
#    try:
#        print(f"Enter the maximum number of repetions (0 <= m <= 6) allowed per colour")
#        m = int(input("m = "))
#        if m > 6:
#            print("Cube has only 6 faces. So, maximum repetion per colour should be less than six\n")
#            continue
#
#        break
#
#    except ValueError:
#        print("Please enter a valid integer")

poly_Z = expand_polynomials(c)
result = count_cube_colorings(poly_Z,m)
print(f"Number of distinct colorings with c={c}, max repetions per color m={m}: {result}")
