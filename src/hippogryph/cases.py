# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from .meshblock import Block, Box, Mesh
from .grid import Uniform, Geometric, Composite

def backward_step(M:int) -> Mesh:
    N = 2*M
    block = Block('domain')
    boxN = Box(ni=17*N, nj=M, block=block, left_label='inflow', right_label='outflow',
               up_label='north')
    boxS = Box(ni=17*N, nj=M, block=block, left_label='south', right_label='outflow',
               down_label='south')
    mesh = Mesh.from_array('BFS', [boxS, boxN], shape=(1,2))

    mesh.index()

    ygrid = Uniform.from_intervals(1.0, mesh.nj, shift=-0.75)
    xunif = Uniform.from_delta(ygrid.delta, 16*N)
    xstretch = Geometric.from_delta(xunif.delta, 16, N)
    xgrid = Composite([xunif, xstretch])

    mesh.mesh(xgrid=xgrid, ygrid=ygrid)

    return mesh

def tee_junction(N:int) -> Mesh:
    
    # parser.add_argument('-o', '--output', dest='output', action='store',
    #                     default='tjunct.exo', help='name of output Exodus II files to be write')
    # parser.add_argument('-N', '--number', dest='N', action='store',
    #                     default=32, help='number of intervals across the channel', type=positive_int)
    # parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
    #                     default=False, help='make lots of noise')
    # parser.add_argument('--no-merge', dest='no_merge', action='store_true',
    #                     default=False, help='keep separate blocks')
    # parser.add_argument('-b', metavar='b', dest='length', action='store',
    #                     default=14, help='sets the branch lengths to L = b*H', type=positive_even_int)
    # parser.add_argument('-i', metavar='i', dest='inlet', action='store',
    #                     default=3, help='sets the inlet length to W = i*H', type=positive_even_int)
    # parser.add_argument('-H', dest='H', action='store', type=positive_float,
    #                     default=1, help='set the height of the inlet, defaults to 1')

    inlet_div_H = 3 # length of the inlet in channel heights
    branch_div_H = 14 # length of the branches in channel heights
    H = 1.0 # height of the inlet channel

    half = int(0.5 * branch_div_H)
    W = inlet_div_H * H

    block = Block('domain')
    inlet = Box(ni=W*N, nj=N, block=block, left_label='inflow', up_label='inlet_north', down_label='south')
    junction = Box(ni=N, nj=N, block=block, down_label='south')
    main0 = Box(ni=half*N, nj=N, block=block, down_label='south', up_label='main_north')
    main1 = Box(ni=N, nj=N, block=block, down_label='south', right_label='east_outflow', up_label='main_north')
    branch0 = Box(ni=N, nj=half*N, block=block, left_label='branch_west', right_label='branch_east')
    branch1 = Box(ni=N, nj=N, block=block, left_label='branch_west', right_label='branch_east', up_label='north_outflow')

    mesh = Mesh.from_array('T-junction', [inlet, junction, main0, main1, None, branch0, None, None, None, branch1, None, None], shape=(4,3))

    # for box in mesh.primitives:
    #     print(box.i, box.j)

    mesh.index()

    delta = H/N
    xunif = Uniform.from_delta(delta, inlet.ni + junction.ni + main0.ni)
    xstretch = Geometric.from_delta(xunif.delta, 7*H, main1.ni)
    xgrid = Composite([xunif, xstretch])
    yunif = Uniform.from_delta(delta, inlet.nj + branch0.nj, shift=-0.5*H)
    ystretch = Geometric.from_delta(yunif.delta, 7*H, branch1.ni)
    ygrid = Composite([yunif, ystretch])

    mesh.mesh(xgrid=xgrid, ygrid=ygrid)

    return mesh
