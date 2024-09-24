# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from .meshblock import Mesh, Box, Block, DimensionalityError
from .grid import Uniform, Geometric, Composite

def half_channel(x:float=1.0, y:float=1.0, z:float=0.0, ni:int=32, nj:int=32, nk:int=0, dy:float|None=None,
                 down_label:str|None=None, up_label:str|None=None, left_label:str|None=None,
                 right_label:str|None=None, front_label:str|None=None, back_label:str|None=None,
                 domain_label:str|None=None) -> Block:
    """
    Generate a half-channel grid
    """
    (up_label, down_label, 
     left_label, right_label, 
     front_label, back_label,
     domain_label) = (el[1] if el[0] is None else el[0] for el in zip([up_label, down_label,
                                                                       left_label, right_label,
                                                                       front_label, back_label,
                                                                       domain_label],
                                                                       ['top_wall', 'bottom_wall',
                                                                        'inflow', 'outflow',
                                                                        'front_wall', 'back_wall',
                                                                        'domain']))
    if nk == 0 and z > 0.0:
        raise DimensionalityError('Channel grid specifies zero z-direction cells but non-zero z length')
    if nk > 0 and  z == 0.0:
        raise DimensionalityError('Channel grid specifies zero z length but non-zero z-direction cells')

    mesh = Mesh('half-channel')
    if nk == 0:
        box = Box(ni=ni, nj=nj, nk=nk, element_set=domain_label, left_label=left_label, right_label=right_label,
                  up_label=up_label, down_label=down_label)
    else:
        box = Box(ni=ni, nj=nj, nk=nk, element_set=domain_label, left_label=left_label, right_label=right_label,
                  up_label=up_label, down_label=down_label, front_label=front_label,
                  back_label=back_label)

    block = mesh.new_block(box)

    block.index()

    xgrid = Uniform.from_intervals(x, block.ni)
    if dy is None:
        ygrid = Uniform.from_intervals(y, block.nj, shift=-0.5*y)
    else:
        ygrid = Geometric.from_delta(dy, y, nj)
    #ygrid = Uniform.from_intervals(y, block.nj, shift=-0.5*y)
    zgrid = None
    if nk != 0:
        zgrid = Uniform.from_intervals(z, block.nk, shift=-0.5*z)

    block.mesh(xgrid=xgrid, ygrid=ygrid, zgrid=zgrid)

    return mesh

def channel(x:float=1.0, y:float=1.0, z:float=0.0, ni:int=32, nj:int=32, nk:int=0,
            down_label:str|None=None, up_label:str|None=None, left_label:str|None=None,
            right_label:str|None=None, front_label:str|None=None, back_label:str|None=None,
            domain_label:str|None=None) -> Block:
    """
    Generate a channel grid
    """
    (up_label, down_label,
     left_label, right_label,
     front_label, back_label,
     domain_label) = (el[1] if el[0] is None else el[0] for el in zip([up_label, down_label,
                                                                       left_label, right_label,
                                                                       front_label, back_label,
                                                                       domain_label],
                                                                      ['top_wall', 'bottom_wall',
                                                                       'inflow', 'outflow',
                                                                       'front_wall', 'back_wall',
                                                                       'domain']))


    if nk == 0 and z > 0.0:
        raise DimensionalityError('Channel grid specifies zero z-direction cells but non-zero z length')
    if nk > 0 and  z == 0.0:
        raise DimensionalityError('Channel grid specifies zero z length but non-zero z-direction cells')

    mesh = Mesh('channel')
    if nk == 0:
        box = Box(ni=ni, nj=nj, nk=nk, element_set=domain_label, left_label=left_label, right_label=right_label,
                  up_label=up_label, down_label=down_label)
    else:
        box = Box(ni=ni, nj=nj, nk=nk, element_set=domain_label, left_label=left_label, right_label=right_label,
                  up_label=up_label, down_label=down_label, front_label=front_label,
                  back_label=back_label)

    block = mesh.new_block(box)

    block.index()

    xgrid = Uniform.from_intervals(x, block.ni)
    ygrid = Uniform.from_intervals(y, block.nj, shift=-0.5*y)
    zgrid = None
    if nk != 0:
        zgrid = Uniform.from_intervals(z, block.nk, shift=-0.5*z)

    block.mesh(xgrid=xgrid, ygrid=ygrid, zgrid=zgrid)

    return mesh

def backward_step(M:int) -> Block:
    """
    Generate a backward-facing step grid
    """
    mesh = Mesh('channel')
    N = 2*M
    boxN = Box(ni=17*N, nj=M, element_set='domain', left_label='inflow', right_label='outflow',
               up_label='north')
    boxS = Box(ni=17*N, nj=M, element_set='domain', left_label='south', right_label='outflow',
               down_label='south')
    block = mesh.block_from_list('block-1', [boxS, boxN], shape=(1,2))

    mesh.index()

    ygrid = Uniform.from_intervals(1.0, block.nj, shift=-0.75)
    xunif = Uniform.from_delta(ygrid.delta, 16*N)
    xstretch = Geometric.from_delta(xunif.delta, 16, N)
    xgrid = Composite([xunif, xstretch])

    block.mesh(xgrid=xgrid, ygrid=ygrid)

    return mesh

def tee_junction(N:int) -> Block:
    """
    Generate a tee-junction grid
    """
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

    mesh = Mesh('T-junction')

    eblock = mesh.new_element_set('domain')
    inlet = Box(ni=inlet_div_H*N, nj=N, element_set=eblock, left_label='inflow',
                up_label='inlet_north', down_label='south')
    junction = Box(ni=N, nj=N, element_set=eblock, down_label='south')
    main0 = Box(ni=half*N, nj=N, element_set=eblock, down_label='south',
                up_label='main_north')
    main1 = Box(ni=N, nj=N, element_set=eblock, down_label='south',
                right_label='east_outflow', up_label='main_north')
    branch0 = Box(ni=N, nj=half*N, element_set=eblock, left_label='branch_west',
                  right_label='branch_east')
    branch1 = Box(ni=N, nj=N, element_set=eblock, left_label='branch_west',
                  right_label='branch_east', up_label='north_outflow')

    block = mesh.block_from_list('block-1', [inlet, junction, main0, main1,
                                             None, branch0, None, None,
                                             None, branch1, None, None], shape=(4,3))

    block.index()

    delta = H/N
    xunif = Uniform.from_delta(delta, inlet.ni + junction.ni + main0.ni)
    xstretch = Geometric.from_delta(xunif.delta, 7*H, main1.ni)
    xgrid = Composite([xunif, xstretch])
    yunif = Uniform.from_delta(delta, inlet.nj + branch0.nj, shift=-0.5*H)
    ystretch = Geometric.from_delta(yunif.delta, 7*H, branch1.ni)
    ygrid = Composite([yunif, ystretch])

    block.mesh(xgrid=xgrid, ygrid=ygrid)

    return mesh
