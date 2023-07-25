# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import hippogryph

def main(inputs = None):
    if inputs is None:
        inputs = {}
    factor = hippogryph.single_sided_vinokur(1.0e-6, 1, 16)
    print(factor)
    for i in range(1, 17):
        print(i, hippogryph.vkruh(factor, 1.0, i, 16))

    print()
    geom = hippogryph.Geometric(1.2, 1.0e-3, 1.0, 16)
    print(geom.s(16))
    delta = 1.0e-3/geom.s(16)
    geom = hippogryph.Geometric(1.2, delta, 1.0, 16)
    print(geom.s(16))
    factor = hippogryph.single_sided_geometric(delta, 16)
    print(factor)
    print(hippogryph.geometric(factor, delta, 16))

    print()
    factor = hippogryph.single_sided_geometric(1.0e-3, 16)
    print(factor)
    for i in range(20):
        print(i, hippogryph.geometric(factor, 1.0e-3, i), 1.0e-3*(1.0 - factor**i)/(1.0 - factor))
    print((factor + 1) * 1.0e-3)

    print()
    factor = hippogryph.single_sided_geometric(1.0e-6, 16, init=200.0)
    print(factor)
    print(hippogryph.geometric(factor, 1.0e-6, 16))
    for i in range(20):
        print(i, hippogryph.geometric(factor, 1.0e-6, i))