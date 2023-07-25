import argparse
import exodusii
import hippogryph as hpg

def positive_even_int(string):
    value = int(string)
    if value < 0 or value % 2 != 0:
        msg = "%r is not an even, positive integer" % string
        raise argparse.ArgumentTypeError(msg)
    return value

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate backward facing step mesh.')
    parser.add_argument('-o', '--output', dest='output', action='store',
                        default='bfs.exo', help='name of output Exodus II files to be write')
    parser.add_argument('-N', '--number', dest='N', action='store',
                        default=32, help='number of intervals across the channel', type=positive_even_int)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='make lots of noise')
    #parser.add_argument('-l', '--list', dest='list', action='store_true',
    #                    default=False, help='list objects and exit')
    #parser.add_argument('-c', '--contains', dest='contains', action='store',
    #                    default='', help='simple filter using "in" for matching')
    #parser.add_argument('-e', '--exactly', dest='exactly', action='append',
    #                    help='simple filter using exact matching')

    args = parser.parse_args()

    N = args.N
    M = int(N/2)

    block = hpg.Block('domain')
    boxN = hpg.Box(ni=17*N, nj=M, block=block, left_label='inflow', right_label='outflow',
                   up_label='north')
    boxS = hpg.Box(ni=17*N, nj=M, block=block, left_label='south', right_label='outflow',
                   down_label='south')
    mesh = hpg.Mesh.from_array('BFS', (1,2), [boxS, boxN])

    for box in mesh.primitives:
        print(box.i, box.j)

    mesh.build()

    xgrid = hpg.Uniform.from_intervals(32.0, mesh.ni)
    ygrid = hpg.Uniform.from_intervals(1.0, mesh.nj)

    mesh.apply(xgrid=xgrid, ygrid=ygrid)

    mesh.save(args.output)