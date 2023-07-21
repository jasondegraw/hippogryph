# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause
import math

def vkruh(factor: float, L: float, i: float, I: int) -> float:
    return L * (1.0 + math.tanh(factor * (i / I - 1.0)) / math.tanh(factor))

def single_sided_factor(ds: float, i: float, I: int, 
                        tolerance: float = 1.0e-14,
                        max_iterations: int = 100,
                        output = print) -> float:
    """
    single_sided_factor - Given a spacing, determine the required stretching factor

    Vinokur's one-sided stretching function between 0 and 1 is

            ds = 1 + tanh(factor * (i / I - 1)) / tanh(factor)

    We apply Newton's Method to find the stretching factor given ds, i, and
    I. The stretching function can be reorganized as

             (ds - 1)tanh(factor) = tanh(factor * (i / I - 1))
                C2 * tanh(factor) = tanh(C1 * factor)

    This only has solutions if |C2| > |C1|.  To see this, note that the
    limit of the LHS is C2, and the limit of the RHS is -1.  Both C1 and C2
    are negative.  The slopes at zero are C2 and C1. respectively, so the
    only possible way that the two curves can intersect is if |C2| > |C1|.
    Approximating tanh(C1*delta) as the line C1*delta for small delta, and
    approximating  C2*tanh(delta) as C2 for large delta, we obtain an first
    guess for the answer as

                               C2 = C1*factor

    Newton's method can now be used to determine the solution.
    """

    C1 = float(i) / float(I) - 1.0
    C2 = ds - 1.0
    
    output("Vinokur h Stretching Factor Solution ---------------+")
    output("  ds = % .8e                              |" % ds)
    output("   I = % .4e                                  |" % I)
    output(" tol = % .3e, itermax = %5d                  |" % (tolerance, max_iterations))
    output("----------------------------------------------------+")
    
    # Solve only for positive ds
    if ds <= 0.0:
       output("No solution for negative ds.                        |")
       output("----------------------------------------------------+")
       return None

    # Only have a nontrivial solution for |C2| > |C1|
    if (C2 > 0.0) or (C2 > C1):
        output("No solution for given inputs.                       |")
        output("----------------------------------------------------+")
        return None
    
    # As an initial guess, take the intersection of the lines
    #      y = 1.0   (limit of tanh(d*C1)
    #      y = C2*d  (approximation of C2*tanh(d))
    #      y = C2    (limit of C2*tanh(d))
    #      y = C1*d  (approximation of tanh(C1*d))

    factor = C1 / C2
    f = math.tanh(factor * C1) - C2 * math.tanh(factor)
    output(" iter          factor                   f           |")
    output("----- ---------------------- ---------------------- |")
    output("%5d % .15e % .15e |" % (1, factor, f))

    for iter in range(2, max_iterations+1):
        s1 = 1.0 / math.cosh(factor * C1)
        s2 = 1.0 / math.cosh(factor)
        fp = C1 * s1 * s1 - C2 * s2 * s2
        factor -= f / fp
        f = math.tanh(factor * C1) - C2 * math.tanh(factor)

        output("%5d % .15e % .15e |" % (iter, factor, f))

        if abs(f) <= tolerance:
            output("----------------------------------------------------+")
            break
    else:
        output("Failed to converge.                                 |")
        output("----------------------------------------------------+")
        return None
    
    return factor

class Uniform:
    def __init__(self, delta: float, L: float, N: int):
        self.delta = delta
        self.L = L
        self.N = N

    def s(self, i:float) -> float:
        return i * self.delta
    
    @classmethod
    def from_delta(cls, delta: float, N: int):
        L = delta * N
        return cls(delta, L, N)
    
    @classmethod
    def from_intervals(cls, L: float, N: int):
        delta = L / float(N)
        return cls(delta, L, N)

class VinokurSingleSided:
    def __init__(self, factor: float, L: float, N: int):
        self.factor = factor
        self.N = N
        self.I = N # The number of intervals is also where the function ends
        self.L = L

    def s(self, i: float) -> float:
        return vkruh(self.factor, self.L, i, self.I)
    
    @classmethod
    def from_delta(cls, delta: float, L: float, i: float, I: int, tolerance: float = 1.0e-14, max_iterations: int = 100,
                   output=print):
        ds = delta / L # Rescale
        factor = single_sided_factor(ds, i, I, tolerance=tolerance, max_iterations=max_iterations, output=output)
        if factor is None:
            return None
        return cls(factor, L, I)
    
class Composite:
    def __init__(self, grids=None):
        self.grids = grids
        if grids is None:
            self.grids = []
        self.L = 0.0
        self.N = 0
        self.intervals = []
        for grid in self.grids:
            self.N += grid.N
            self.L += grid.L
            self.intervals.append(self.N)

    def s(self, i: float) -> float:
        for k, N in enumerate(self.intervals):
            if i <= N:
                return self.grids[k].s(i) # Not correct
        return None