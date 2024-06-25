from time import perf_counter

import sys

sys.path.insert(1, "/Users/mpham26/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/laughing-umbrella/xx_synthesis/monodromy")


import monodromy

from monodromy.coordinates import monodromy_alcove, monodromy_alcove_c2, monodromy_to_positive_canonical_polytope, rho_reflect
from monodromy.elimination import cylinderize, project
from monodromy.polytopes import ConvexPolytope, Polytope
from monodromy.static import qlr_polytope

from itertools import count

biswas_relations = (qlr_polytope
    # enlarge to the pu_4 version of the QLR relations
    .union(rho_reflect(qlr_polytope, [0, 7, 8, 9]))
    # constrain in- and out-coordinates to the appropriate alcove
    .intersect(cylinderize(monodromy_alcove, [0, 1, 2, 3], 10))
    .intersect(cylinderize(monodromy_alcove_c2, [0, 7, 8, 9], 10))
)

# constrain interaction coordinates to be of XY-type
biswas_relations = biswas_relations.intersect(Polytope(convex_subpolytopes=[
    ConvexPolytope(
        inequalities=[[1, 0, 0, 0, -4, 0, 0, 0, 0, 0]],
        equalities=[
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # x2 == 0
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # x3 == 0
        ]
    )
]))

# switch to canonical coordinates
biswas_relations = monodromy_to_positive_canonical_polytope(
    biswas_relations, coordinates=[0, 1, 2, 3])
biswas_relations = monodromy_to_positive_canonical_polytope(
    biswas_relations, coordinates=[0, 4, 5, 6])
biswas_relations = monodromy_to_positive_canonical_polytope(
    biswas_relations, coordinates=[0, 7, 8, 9])

# reduce the biswas relations to have following coordinates:
# k a1 a2 a3 beta b1 b2 b3
biswas_relations = biswas_relations.reduce()
biswas_relations = project(biswas_relations, 6).reduce()
biswas_relations = project(biswas_relations, 5).reduce()
    
xy_polytope = monodromy.static.examples.identity_polytope
polytope_list = []
for n in range(4):
    
    print(f"Working on an XY interaction sequence of length {n}...")
    
    # inflate xy_polytope from [*a_coords, *interaction_coords] to [*a_coords, *b_coords, *interaction_coords, beta]
    xy_polytope = cylinderize(
        xy_polytope,
        coordinate_map=[0, 1, 2, 3] + list(range(7, 7 + (n - 1))),
        parent_dimension=1 + 3 + 3 + n,
    ).intersect(cylinderize(
        biswas_relations,
        coordinate_map=[0, 1, 2, 3, -1, 4, 5, 6],
        parent_dimension=1 + 3 + 3 + n,
    ))
    
    # project away the old a-coordinates
    start_time = perf_counter()
    print("Working on the reduction 1/3...", end="")
    xy_polytope = project(xy_polytope, 3).reduce()
    print(f" done.  Took {perf_counter() - start_time} seconds.")

    start_time = perf_counter()
    print("Working on the reduction 2/3...", end="")
    xy_polytope = project(xy_polytope, 2).reduce()
    print(f" done.  Took {perf_counter() - start_time} seconds.")

    start_time = perf_counter()
    print("Working on the reduction 3/3...", end="")
    xy_polytope = project(xy_polytope, 1).reduce()
    print(f" done.  Took {perf_counter() - start_time} seconds.")
    
    # now the old c-coordinates are sitting where the a-coordinates were!
    print("The first three coordinates are the canonical coordinates CAN(x1, x2, x3).")
    print("The remaining coordinates x4, ..., xk are the XY interaction strengths.")
    print(xy_polytope)
    
    polytope_list.append(xy_polytope)

def print_table(data, no_extra_rows, name):
    # Define the headers
    headers = ["1", "a_1", "a_2", "a_3"] + [f"𝛼_{i+1}" for i in range(no_extra_rows)]
    
    data = sorted(data)
    
    # Ensure that the data length matches the number of headers
    if not all(len(row) == len(headers) for row in data):
        raise ValueError("All rows in the data must have the same length as the headers.")
    
    # Calculate the maximum width for each column
    col_widths = [max(len(str(item)) for item in [header] + [row[i] for row in data]) for i, header in enumerate(headers)]
    
    # Define a function to print a single row with given data and column widths
    def print_row(row):
        print(" | ".join(f"{str(item).ljust(col_widths[i])}" for i, item in enumerate(row)))
    
    # Print name
    print(f"{name} ({len(data)} inequalities)\n")
    
    # Print the headers
    print_row(headers)
    print("-+-".join("-" * width for width in col_widths))
    
    # Print each row of data
    for row in data:
        print_row(row)
        
    print("\n")

def classify(ineq_system):
    
    interaction_ineq, CAN_ineq, total_strength_ineq = [], [], []
    slant_ineq, height_ineq, other_ineq = [], [], []
    spade_ineq, club_ineq, diamond_ineq, heart_ineq = [], [], [], []
    for ineq in ineq_system:
        if ineq[1:4] == [0, 0, 0]:
            interaction_ineq.append(ineq)
        elif ineq[4:] == [0 for i in range(len(ineq[4:]))]:
            CAN_ineq.append(ineq)
        elif ineq[1:4] == [-1, -1, -1]:
            total_strength_ineq.append(ineq)
        elif ineq[1:4] == [1, -1, -1]:
            slant_ineq.append(ineq)
        elif ineq[1:4] == [0, 0, -1]:
            height_ineq.append(ineq)
        elif ineq[1:4] == [-1, 0, 0]:
            spade_ineq.append(ineq)
        elif ineq[1:4] == [-1, 1, -1]:
            club_ineq.append(ineq)
        elif ineq[1:4] == [0, 1, 0]:
            diamond_ineq.append(ineq)
        elif ineq[1:4] == [1, 1, -1]:
            heart_ineq.append(ineq)
        else:
            other_ineq.append(ineq)
    
    n = len(ineq)-4
    print_table(interaction_ineq, n, "Interaction Ineq")
    print_table(CAN_ineq, n, "CAN Ineq")
    print_table(total_strength_ineq, n, "Total Strength Ineq")
    print_table(slant_ineq, n, "Slant Ineq")
    print_table(height_ineq, n, "Height Ineq")
    print_table(spade_ineq, n, "Spade Ineq")
    print_table(club_ineq, n, "Club Ineq")
    print_table(diamond_ineq, n, "Diamond Ineq")
    print_table(heart_ineq, n, "Heart Ineq")
    print_table(other_ineq, n, "Other Ineq")
                
    return interaction_ineq, CAN_ineq, total_strength_ineq, slant_ineq, height_ineq, other_ineq

polytope = polytope_list[3].convex_subpolytopes[0].inequalities
interaction_ineq, CAN_ineq, total_strength_ineq, slant_ineq, height_ineq, other_ineq = classify(polytope)