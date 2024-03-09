# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import exodusii
import hippogryph as hpg

def positive_even_int(string):
    value = int(string)
    if value < 0 or value % 2 != 0:
        msg = "%r is not an even, positive integer" % string
        raise argparse.ArgumentTypeError(msg)
    return value

def positive_int(string):
    value = int(string)
    if value < 0:
        msg = "%r is not a positive integer" % string
        raise argparse.ArgumentTypeError(msg)
    return value

def positive_float(string):
    value = float(string)
    if value < 0:
        msg = "%r is not a positive floating point number" % string
        raise argparse.ArgumentTypeError(msg)
    return value

txt = '''Generate T-junction mesh:

             y2 +-------+
                |       |  * H = y1 - y0 = x2 - x1
                |       |  * L = x3 - x2 = y2 - y1
                |       |  * All H length segments are 
                |       |    divided into N equal intervals
                |       |  * d = H/N
                | Block |  * W = x1 - x0 is an integer
                |   4   |    multiple of H
                |       |  * Before/below x2 + L/2 and
                |       |    y1 + L/2, the elements are
                |       |    d x d squares
                |       |  * After/above x2 + L/2 and
                |       |    y1 + L/2, N intervals with
                |       |    d0 = d stretching are used
   y1 +---------+ - - - +-----------------------------------+
      | Block   . Block .           Block                   |
      |   1     .   2   .             3                     |
   y0 +---------+-------+-----------------------------------+
     x0        x1     x2                                  x3   


'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=txt)
    parser.add_argument('-o', '--output', dest='output', action='store',
                        default='tjunct.exo', help='name of output Exodus II files to be write')
    parser.add_argument('-N', '--number', dest='N', action='store',
                        default=32, help='number of intervals across the channel', type=positive_int)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='make lots of noise')
    parser.add_argument('--no-merge', dest='no_merge', action='store_true',
                        default=False, help='keep separate blocks')
    parser.add_argument('-b', metavar='b', dest='length', action='store',
                        default=14, help='sets the branch lengths to L = b*H', type=positive_even_int)
    parser.add_argument('-i', metavar='i', dest='inlet', action='store',
                        default=3, help='sets the inlet length to W = i*H', type=positive_even_int)
    parser.add_argument('-H', dest='H', action='store', type=positive_float,
                        default=1, help='set the height of the inlet, defaults to 1')

    args = parser.parse_args()

    H = args.H
    N = args.N
    half = int(0.5 * args.length)
    Lo2 = half * H
    W = args.inlet * H

    block = hpg.Block('domain')
    inlet = hpg.Box(ni=W*N, nj=N, block=block, left_label='inflow', up_label='inlet_north', down_label='south')
    junction = hpg.Box(ni=N, nj=N, block=block, down_label='south')
    main0 = hpg.Box(ni=half*N, nj=N, block=block, down_label='south', up_label='main_north')
    main1 = hpg.Box(ni=N, nj=N, block=block, down_label='south', right_label='east_outflow', up_label='main_north')
    branch0 = hpg.Box(ni=N, nj=half*N, block=block, left_label='branch_west', right_label='branch_east')
    branch1 = hpg.Box(ni=N, nj=N, block=block, left_label='branch_west', right_label='branch_east', up_label='north_outflow')

    mesh = hpg.Mesh.from_array('T-junction', [inlet, junction, main0, main1, None, branch0, None, None, None, branch1, None, None], shape=(4,3))

    for box in mesh.primitives:
        print(box.i, box.j)

    mesh.build()

    delta = H/N
    xunif = hpg.Uniform.from_delta(delta, inlet.ni + junction.ni + main0.ni)
    xstretch = hpg.Geometric.from_delta(xunif.delta, 7*H, main1.ni)
    xgrid = hpg.Composite([xunif, xstretch])
    yunif = hpg.Uniform.from_delta(delta, inlet.nj + branch0.nj, shift=-0.5*H)
    ystretch = hpg.Geometric.from_delta(yunif.delta, 7*H, branch1.ni)
    ygrid = hpg.Composite([yunif, ystretch])

    mesh.apply(xgrid=xgrid, ygrid=ygrid)

    mesh.save(args.output)

    print(mesh.cell_count)