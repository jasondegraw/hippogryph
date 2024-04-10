# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import numpy
#import exodusii

import hippogryph as hpg

def backward_step(M:int):
    N = 2*M
    block = hpg.Block('domain')
    boxN = hpg.Box(ni=17*N, nj=M, block=block, left_label='inflow', right_label='outflow',
                   up_label='north')
    boxS = hpg.Box(ni=17*N, nj=M, block=block, left_label='south', right_label='outflow',
                   down_label='south')
    mesh = hpg.Mesh.from_array('BFS', [boxS, boxN], shape=(1,2))

    mesh.build()

    ygrid = hpg.Uniform.from_intervals(1.0, mesh.nj, shift=-0.75)
    xunif = hpg.Uniform.from_delta(ygrid.delta, 16*N)
    xstretch = hpg.Geometric.from_delta(xunif.delta, 16, N)
    xgrid = hpg.Composite([xunif, xstretch])

    mesh.apply(xgrid=xgrid, ygrid=ygrid)

    return mesh
