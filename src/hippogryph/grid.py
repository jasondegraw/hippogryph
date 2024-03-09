# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import math
import numbers

def vkruh(factor: float, L: float, i: float, I: int) -> float:
    return L * (1.0 + math.tanh(factor * (i / I - 1.0)) / math.tanh(factor))

def geometric_sum(factor: float, delta: float, i: float) -> float:
    if i == 0:
        return 0.0
    if i == 1:
        return delta
    sum = delta
    for k in range(2, i+1):
        delta *= factor
        sum += delta
    return sum

def geometric(factor: float, delta: float, i: float) -> float:
    if i == 0:
        return 0.0
    if factor == 1:
        return i*delta
    return delta * (1.0 - factor**i)/(1.0 - factor)

def single_sided_geometric(delta: float,
                           I: int, 
                           tolerance: float = 1.0e-14,
                           max_iterations: int = 100,
                           output = print, init=2.0) -> float:
    """
    single_sided_geometric - Given a spacing, determine the required stretching factor

    The geometric one-sided stretching function between 0 and 1 is

            s = sum_k=0,i-1 (delta * factor^k)

    We apply Newton's Method to find the stretching factor given delta and I.The
    stretching function can be reorganized as

            f(factor) = 1.0 + delta - sum_k=1,I-1 (delta * factor^k)

    Then
            f'(factor) = -sum_k=1,I-1 (k * delta * factor^(k-1))
                       = -delta - sum_k=2,I-1 (k * delta * factor^(k-1))
                       = -delta - sum_k=1,I-2 ((k+1) * delta * factor^k)

    There are probably shortcuts, further work is needed. Newton's method can now
    be used to determine the solution.
    """
    output("Geometric Stretching Factor Solution ---------------+")
    output("  ds = % .8e                              |" % delta)
    output("   I = % .4e                                  |" % I)
    output(" tol = % .3e, itermax = %5d                  |" % (tolerance, max_iterations))
    output("----------------------------------------------------+")
    
    # Solve only for positive delta
    if delta <= 0.0:
       output("No solution for negative ds.                        |")
       output("----------------------------------------------------+")
       return None

    factor = init
    if delta > 1.0:
        output("No solution for negative delta > 1.0.               |")
        output("----------------------------------------------------+")
        return None
    #elif I * delta > 1.0:
        #factor /= I * delta
        #f = 1.0
        #factor_power = 1.0
        #for k in range(1,I):
        #    factor_power *= factor
        #    f -= delta * factor_power
    #else:
        #f = 1.0 - (I-1)*delta

    f = delta * (1.0 - factor**I) - 1.0 + factor

    output(" iter          factor                   f           |")
    output("----- ---------------------- ---------------------- |")
    output("%5d % .15e % .15e |" % (1, factor, f))

    for iter in range(2, max_iterations+1):
        #fp = -delta
        #factor_power = 1.0
        #for k in range(1, I-1):
        #    factor_power *= factor
        #    fp -= (k + 1) * delta * factor_power
        #print(fp)

        fp = -delta * I * factor**(I-1) + 1.0

        factor -= f / fp
        #f = 1.0
        #factor_power = 1.0
        #for k in range(1,I):
        #    factor_power *= factor
        #    f -= delta * factor_power

        f = delta * (1.0 - factor**I) - 1.0 + factor

        output("%5d % .15e % .15e |" % (iter, factor, f))

        if abs(f) <= tolerance:
            output("----------------------------------------------------+")
            break
    else:
        output("Failed to converge.                                 |")
        output("----------------------------------------------------+")
        return None
    
    return factor
    

def single_sided_vinokur(ds: float, i: float, I: int, 
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
    def __init__(self, delta: float, L: float, N: int, shift: float = 0.0):
        self.delta = delta
        self.L = L
        self.N = N
        self.shift = shift

    def s(self, i:float) -> float:
        return self.shift + i * self.delta
    
    @classmethod
    def from_json(cls, **kwargs):
        shift = kwargs.get('shift', 0.0)
        N = kwargs.get('N', None)
        try:
            N = int(N)
        except ValueError:
            return None
        delta = kwargs.get('delta', None)
        if isinstance(delta, numbers.Number):
            L = kwargs.get('L', None)
            if isinstance(L, numbers.Number):
                return cls(delta, L, N, shift=shift)
            else:
                return cls.from_delta(delta, N, shift=shift)
        else:
            L = kwargs.get('L', None)
            if isinstance(L, numbers.Number):
                return cls.from_intervals(L, N, shift=shift)
        return None
    
    @classmethod
    def from_delta(cls, delta: float, N: int, shift: float = 0.0):
        L = delta * N
        return cls(delta, L, N, shift=shift)
    
    @classmethod
    def from_intervals(cls, L: float, N: int, shift: float = 0.0):
        delta = L / float(N)
        return cls(delta, L, N, shift=shift)

class VinokurSingleSided:
    def __init__(self, factor: float, L: float, N: int, shift: float = 0.0):
        self.factor = factor
        self.N = N
        self.I = N # The number of intervals is also where the function ends
        self.L = L
        self.shift = shift

    def s(self, i: float) -> float:
        return self.shift + vkruh(self.factor, self.L, i, self.I)
    
    @classmethod
    def from_json(cls, **kwargs):
        shift = kwargs.get('shift', 0.0)
        factor = kwargs.get('factor', None)
        if factor is None: # Should really check that it's a number
            if all(name in kwargs for name in ('delta', 'i', 'N', 'L')):
                tolerance = kwargs.get('tolerance', 1.0e-14)
                max_iterations = kwargs.get('max_iterations', 100)
                return cls.from_delta(kwargs['delta'], kwargs['L'], kwargs['i'], kwargs['N'],
                                       tolerance=tolerance, max_iterations=max_iterations, shift=shift)
        elif all(name in kwargs for name in ('N', 'L')):
            return cls(factor, kwargs['L'], kwargs['N'], shift=shift)
        return None
    
    @classmethod
    def from_delta(cls, delta: float, L: float, i: float, I: int, tolerance: float = 1.0e-14, max_iterations: int = 100,
                   output=print, shift: float = 0.0):
        ds = delta / L # Rescale
        factor = single_sided_vinokur(ds, i, I, tolerance=tolerance, max_iterations=max_iterations, output=output)
        if factor is None:
            return None
        return cls(factor, L, I, shift=shift)
    
class Geometric:
    def __init__(self, factor: float, delta: float, L: float, N: int, shift: float = 0.0):
        self.factor = factor
        self.delta = delta
        self.L = L
        self.N = N
        self.shift = shift

    def s(self, i: float) -> float:
        return self.shift + geometric(self.factor, self.delta, i)
    
    @classmethod
    def from_json(cls, **kwargs):
        shift = kwargs.get('shift', 0.0)
        if all(name in kwargs for name in ('delta', 'N')):
            if all(name in kwargs for name in ('factor', 'L')):
                return cls(kwargs['factor'], kwargs['delta'], kwargs['L'], kwargs['N'], shift=shift)
            L = kwargs.get('L', 1.0)
            tolerance = kwargs.get('tolerance', 1.0e-14)
            max_iterations = kwargs.get('max_iterations', 100)
            return cls.from_delta(kwargs['delta'], L, kwargs['N'], tolerance=tolerance, max_iterations=max_iterations,
                                  shift=shift)
    
    @classmethod
    def from_delta(cls, delta: float, L: float, I: int, tolerance: float = 1.0e-14, max_iterations: int = 100,
                   output=print, shift: float = 0.0):
        ds = delta / L # Rescale
        factor = single_sided_geometric(ds, I, tolerance=tolerance, max_iterations=max_iterations, output=output)
        if factor is None:
            return None
        return cls(factor, delta, L, I, shift=shift)
    
class Composite:
    def __init__(self, grids=None):
        self.grids = grids
        if not grids:
            self.grids = []
        self.L = 0.0
        self.N = 0
        self.intervals = []
        self.shift = grids[0].shift
        for grid in self.grids:
            self.N += grid.N
            grid.shift = self.L
            self.L += grid.L
            self.intervals.append(self.N)

    def s(self, i: float) -> float:
        last = 0
        for k, n in enumerate(self.intervals[:-1]):
            if i <= n:
                return self.shift + self.grids[k].s(i-last)
            last = n
        if i == self.N:
            return self.shift + self.L
        return self.shift + self.grids[-1].s(i-last)