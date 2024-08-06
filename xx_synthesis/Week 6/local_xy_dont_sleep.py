import numpy as np
from datetime import datetime
import sys

sys.path.insert(1, "/Users/minhpham/Documents/Research/laughing-umbrella/xx_synthesis/monodromy")

import monodromy
import numpy as np

from monodromy.coordinates import monodromy_alcove, monodromy_alcove_c2, \
    monodromy_to_positive_canonical_polytope
from monodromy.coverage import rho_reflect
from monodromy.elimination import cylinderize, project
from monodromy.io.base import ConvexPolytopeData, PolytopeData
from monodromy.polytopes import ConvexPolytope, make_convex_polytope, Polytope

from monodromy.static.examples import empty_polytope
from monodromy.static.qlr_table import qlr_polytope

# STRENGTH POLYTOPE
strength_polytope = make_convex_polytope([
    # k s+ s1 s2 beta
    [999, -1,  0,  0,  0],  # 999 - s+ ≥ 0
    [  0,  1, -1, -1,  0],  # s+ - s1 - s2 ≥ 0
    [  0,  0,  1, -1,  0],  # s1 - s2 ≥ 0
    [  0,  0,  0,  1,  0],  # s2 ≥ 0
    [  1,  0, -3,  0,  0],  # 1/3 - s1 ≥ 0
    [  1,  0,  0,  0, -3],  # 1/3 - beta ≥ 0
    [  0,  0,  0,  0,  1],  # beta ≥ 0
], name="I")

# A POLYTOPE
a_polytope = make_convex_polytope([
    # k ah al af s+ s1 s2
    [0,  1, -1,  0, 0, 0, 0],  # ah - al ≥ 0
    [0,  0,  1,  0, 0, 0, 0],  # al ≥ 0
    [0,  0,  0,  1, 0, 0, 0],  # af ≥ 0
    [1, -1, -1,  0, 0, 0, 0],  # 1 - ah - al ≥ 0
    [1, -1,  0, -1, 0, 0, 0],  # 1 - ah - af ≥ 0
], name="A alcove").intersect(
    # MGC inequalities
    empty_polytope.union(
        # unreflected
        make_convex_polytope([
            # strength bound
            [0, -1, -1, -1, 2, 0, 0],  # - ah - al - af + 2 s+ ≥ 0
        ], name="A unreflected").intersect(
            # slant bound
            empty_polytope.union(make_convex_polytope([
                [0,  1, -1, -1, 2, -2, -2]  # ah - al - af + 2s+ - 2 s1 - 2s2 ≥ 0
            ], name="ah slant")).union(make_convex_polytope([
                [0, -1, -1,  1, 2, -2, -2]  # - ah - al + af + 2s+ - 2 s1 - 2s2 ≥ 0
            ], name="af slant"))
        ).intersect(
            # spade bound
            empty_polytope.union(make_convex_polytope([
                [0,  -1, 0, 0, 1, 0, 0]  # -ah + s+ ≥ 0
            ], name="ah spade")).union(make_convex_polytope([
                [0, 0, 0,  -1, 1, 0, 0]  # -af + s+ ≥ 0
            ], name="af spade"))
        ).intersect(
            # club bound
            empty_polytope.union(make_convex_polytope([
                [0,  1, -1, -1, 2, -2, 0]  # ah - al - af + 2s+ - 2s1 ≥ 0
            ], name="ah club")).union(make_convex_polytope([
                [0, -1, 1,  -1, 2, -2, 0]  # -ah + al -af + 2s+ - 2s1 ≥ 0
            ], name="al club")).union(make_convex_polytope([
                [0, -1, -1,  1, 2, -2, 0]  # -ah - al + af + 2s+ - 2s1 ≥ 0
            ], name="af club"))
        ).intersect(
            # heart bound
            empty_polytope.union(make_convex_polytope([
                [0,  1, -1, 1, 2, -4, 0]  # ah - al + af + 2s+ - 4s1 ≥ 0
            ], name="al heart")).union(make_convex_polytope([
                [0, 1, 1,  -1, 2, -4, 0]  # ah + al - af + 2s+ - 4s1 ≥ 0
            ], name="af heart"))
        )
    ).union(
        # reflected
        make_convex_polytope([
            # slant bound
            [1, -1, -1, -1, 2, -2, -2],  # 1 - ah - al - af + 2s+ - 2s1 - 2s2 ≥ 0
        ], name="A reflected").intersect(
            # strength bound
            empty_polytope.union(make_convex_polytope([
                [0,  1, -1, -1, 2, -4, 0]  # ah - al - af + 2s+ - 4s1 ≥ 0
            ], name="ah strength")).union(make_convex_polytope([
                [0, -1, -1,  1, 2, -4, 0]  # - ah - al + af + 2s+ - 4s1 ≥ 0
            ], name="af slant"))
        ).intersect(
            # spade bound
            empty_polytope.union(make_convex_polytope([
                [-1,  1, 0, 0, 1, 0, 0]  # -1 + ah + s+ ≥ 0
            ], name="ah spade")).union(make_convex_polytope([
                [-1, 0, 0,  1, 1, 0, 0]  # -1 + af + s+ ≥ 0
            ], name="af spade"))
        ).intersect(
            # club bound
            empty_polytope.union(make_convex_polytope([
                [-1, 1, -1,  1, 2, -2, 0]  # -1 + ah - al + af + 2s+ - 2s1 ≥ 0
            ], name="al club")).union(make_convex_polytope([
                [-1,  1, 1, -1, 2, -2, 0]  # -1 + ah + al - af + 2s+ - 2s1 ≥ 0
            ], name="af club"))
        ).intersect(
            # heart bound
            empty_polytope.union(make_convex_polytope([
                [1,  1, -1, -1, 2, -4, 0]  # 1 + ah - al - af + 2s+ - 4s1 ≥ 0
            ], name="ah heart")).union(make_convex_polytope([
                [1, -1, 1,  -1, 2, -4, 0]  # 1 - ah + al - af + 2s+ - 4s1 ≥ 0
            ], name="al heart")).union(make_convex_polytope([
                [1, -1, -1,  1, 2, -4, 0]  # 1 - ah - al + af + 2s+ - 4s1 ≥ 0
            ], name="af heart"))
        )
    ).intersect(
        # frustrum bound
        empty_polytope.union(make_convex_polytope([
            [0, 0, -1,  0, 1, -1, 0]  # -al + s+ - s1 ≥ 0
        ], name="al frustrum")).union(make_convex_polytope([
            [0, 0,  0, -1, 1, -1, 0]  # -af + s+ - s1 ≥ 0
        ], name="af frustrum"))
    ).intersect(
        # diamond bound
        empty_polytope.union(make_convex_polytope([
            [0, 1, 0,  0, 1, -2, 0]  # ah + s+ - 2s1 ≥ 0
        ], name="ah diamond")).union(make_convex_polytope([
            [0, 0,  0, 1, 1, -2, -0]  # af + s+ - 2s1 ≥ 0
        ], name="af diamond"))
    )
)

# B POLYTOPE
b_polytope = make_convex_polytope([
    # k b1 b2 b3 s+ s1 s2 beta
    [0,  1, -1,  0, 0, 0, 0, 0],  # b1 - b2 ≥ 0
    [0,  0,  1, -1, 0, 0, 0, 0],  # b2 - b3 ≥ 0
    [0,  0,  0,  1, 0, 0, 0, 0],  # b3 ≥ 0
    [1, -1, -1,  0, 0, 0, 0, 0],  # 1 - b1 - b2 ≥ 0
], name="B alcove").intersect(
    # MGC inequalities
    empty_polytope.union(make_convex_polytope([
        # strength
        [0, -1, -1, -1, 2,  0,  0,  2],  # - b1 - b2 - b3 + 2s+ + 2beta ≥ 0
        # slant
        [0,  1, -1, -1, 2, -2,  -2,  2],  # b1 - b2 - b3 + 2s+ - 2s1 - 2s2 + 2beta ≥ 0
        [0,  1, -1, -1, 2,  -2,  0, 0],  # b1 - b2 - b3 + 2s+ - 2s1 ≥ 0
        # frustrum
        [0,  0,  0, -1, 1, -1, 0,  1],  # - b3 + s+ - s1 + beta ≥ 0
        [0,  0,  0, -1, 1, 0,  0,  0],  # - b3 + s+ ≥ 0
        # spade
        [0, -1, 0, 0, 1, 0, 0, 1], # - b1 + s+ + beta ≥ 0
        # club
        [0, -1, 1, -1, 2, -2, 0, 2], # - b1 + b2 - b3 + 2s+ -2s1 + 2beta ≥ 0
        [0, -1, 1, -1, 2, 0, 0, 0], # - b1 + b2 - b3 + 2s+ ≥ 0
        # diamond
        [0, 0, 1, 0, 1, -2, 0, 1], # b2 + s+ - 2s1 + beta ≥ 0
        [0, 0, 1, 0, 1, 0, 0, -1], # b2 + s+ - beta ≥ 0
        # heart
        [0, 1, 1, -1, 2, -4, 0, 2], # b1 + b2 - b3 + 2s+ - 4s1 + 2beta ≥ 0
        [0, 1, 1, -1, 2, 0, 0, -2], # b1 + b2 - b3 + 2s+ - 2beta ≥ 0
    ], name="B unreflected")).union(make_convex_polytope([
        # strength
        [-1 , 1, -1, -1, 2,  0,  0,  2],  # -1 + b1 - b2 - b3 + 2s+ + 2beta ≥ 0
        # slant
        [ 1, -1, -1, -1, 2,  -2,  -2, 2],  # 1 - b1 - b2 - b3 + 2s+ - 2s1 - 2s2 + 2beta ≥ 0
        [ 1, -1, -1, -1, 2, -2,  0,  0],  # 1 - b1 - b2 - b3 + 2s+ - 2s1 ≥ 0
        # frustrum
        [ 0,  0,  0, -1, 1, -1, 0,  1],  # - b3 + s+ - s1 + beta ≥ 0
        [ 0,  0,  0, -1, 1, 0,  0,  0],  # - b3 + s+ ≥ 0
        # spade
        [-1, 1, 0, 0, 1, 0, 0, 1],  # -1 + b1 + s+ + beta ≥ 0
        # club
        [-1, 1, 1, -1, 2, -2, 0, 2], # -1 + b1 + b2 - b3 + 2s+ - 2s1 + 2beta ≥ 0
        [-1, 1, 1, -1, 2, 0, 0, 0], # -1 + b1 + b2 - b3 + 2s+ ≥ 0
        # diamond
        [0, 0, 1, 0, 1, -2, 0, 1], # b2 + s+ - 2s1 + beta ≥ 0
        [0, 0, 1, 0, 1, 0, 0, -1], # b2 + s+ - beta ≥ 0
        # heart
        [1, -1, 1, -1, 2, -4, 0, 2], # 1 - b1 + b2 - b3 + 2s+ - 4s1 + 2beta ≥ 0
        [1, -1, 1, -1, 2, 0, 0, -2], # 1 - b1 + b2 - b3 + 2s+ - 2beta ≥ 0
    ], name="B reflected"))
)

# INTERFERENCE POLYTOPE
interference_polytope = Polytope(convex_subpolytopes=[
    # af = b1
    ConvexPolytope([
        # k ah al af b1 b2 b3 beta
        [0, -1, -1, 0,  0,  1,  1,  2],  # - ah - al + b2 + b3 + 2beta ≥ 0
        [2, -1, -1, 0,  0, -1, -1, -2],  # 2 - ah - al - b2 - b3 - 2beta ≥ 0
        [0,  1,  1, 0,  0, -1, -1,  2],  # ah + al - b2 - b3 + 2beta ≥ 0
    ], equalities=[
        [0,  0,  0, 1, -1,  0,  0,  0],  # af = b1
        [0, 1, -1, 0, 0, -1, 1, 0],  # ah - al = b2 - b3
    ], name="AF=B1"),
    # OR af = b2
    ConvexPolytope([
        # k ah al af b1 b2 b3 beta
        [0, -1, -1, 0,  1,  0,  1,  2],  # - ah - al + b1 + b3 + beta ≥ 0
        [2, -1, -1, 0, -1,  0, -1, -2],  # 2 - ah - al - b1 - b3 - beta ≥ 0
        [0,  1,  1, 0, -1,  0, -1,  2],  # ah + al - b1 - b3 + beta ≥ 0
    ], equalities=[
        [0,  0,  0, 1,  0, -1,  0,  0],  # af = b2
        [0, 1, -1, 0, -1, 0, 1, 0],  # ah - al = b1 - b3
    ], name="AF=B2"),
    # OR af = b3
    ConvexPolytope([
        # k ah al af b1 b2 b3 beta
        [0, -1, -1, 0,  1,  1,  0,  2],  # - ah - al + b1 + b2 + beta ≥ 0
        [2, -1, -1, 0, -1, -1,  0, -2],  # 2 - ah - al - b1 - b2 - beta ≥ 0
        [0,  1,  1, 0, -1, -1,  0,  2],  # ah + al - b1 - b2 + beta ≥ 0
    ], equalities=[
        [0,  0,  0, 1,  0,  0, -1,  0],  # af = b3
        [0, 1, -1, 0, -1, 1, 0, 0],  # ah - al = b1 - b2
    ], name="AF=B3"),
])

def regenerate_xy_solution_polytopes():
    """
    Recalculates the partition of b-coordinates into regions with solvable lifts
    with the interference relations.

    Recreates the pair (xy_region_polytope, xy_lift_polytope).

    NOTE: This routine amounts to a computer-calculated _proof_ of the main
          local theorem, which it checks as an assertion.
    """
    
    print("Start time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Checking main local theorem.")

    raw_polytope = (
        # k s+ s1 s2 beta
        cylinderize(strength_polytope, [0, 7, 8, 9, 10], 11)
        # k ah al af s+ s1 s2
        .intersect(cylinderize(a_polytope, [0, 1, 2, 3, 7, 8, 9], 11))
        # k b1 b2 b3 s+ s1 s2 beta
        .intersect(cylinderize(b_polytope, [0, 4, 5, 6, 7, 8, 9, 10], 11))
        # k ah al af b1 b2 b3 beta
        .intersect(cylinderize(interference_polytope, [0, 1, 2, 3, 4, 5, 6, 10], 11))
    )

    constrained_polytope = raw_polytope.reduce()
    # project away the a polytope: af, then al, then ah.
    constrained_polytope = project(constrained_polytope, 3).reduce()
    constrained_polytope = project(constrained_polytope, 2).reduce()
    constrained_polytope = project(constrained_polytope, 1).reduce()

    # compare with the original b polytope
    big_polytope = b_polytope.intersect( 
        cylinderize(strength_polytope, [0, 4, 5, 6, 7], 8)
    ).reduce()
    assert constrained_polytope.contains(big_polytope), "Not true"
    
    print("Done.")
    print("End time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    return constrained_polytope

poly = [regenerate_xy_solution_polytopes()]

print(poly)