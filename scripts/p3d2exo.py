# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import hippogryph as hpg
import os
import struct
import numpy as np

def devnull(mesg):
    pass

def node_index(i, j, k, ni, nj, nk):
    return i + 1 + j * ni + k * ni * nj

def read_plotXd(input_file:str, status=devnull):
    # Just assume it exists
    status('Attempting to read "%s" as a 3D, multiblock, unformatted file...' % input_file)
    try:
        return hpg.plot3d.read_plot3D(input_file)
    except struct.error:
        # Nope, not 3D, multiblock, unformatted
        pass
    status('Attempting to read "%s" as a 3D, multiblock, formatted file...' % input_file)
    try:
        return hpg.plot3d.read_plot3D(input_file, binary=False)
    except IndexError:
        # Not 3D, multiblock, formatted
        pass
    status('Attempting to read "%s" as a 2D, multiblock, unformatted file...' % input_file)
    try:
        return hpg.plot3d.read_plot2D(input_file)
    except struct.error:
        # Not 2D, multiblock, unformatted
        pass
    return hpg.plot3d.read_plot2D(input_file, binary=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate backward facing step mesh.')
    parser.add_argument('-o', '--output', dest='output', action='store',
                        default='p3d2exo.exo', help='name of output Exodus II files to be write')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=False, help='make lots of noise')
    parser.add_argument('P3D', action='store', help='Input plot3d file to process')

    args = parser.parse_args()

    if not os.path.exists(args.P3D):
        exit('Unable to find input file "%s"' % args.P3D)

    input_file = args.P3D
    output_file = args.output
    verbose = args.verbose

    status = devnull
    if verbose:
        status = print

    two_dimensional = False

    blocks = read_plotXd(input_file, status=status)

    status('Read %d block(s)' % len(blocks))
    if verbose:
        for block in blocks:
            nk = max(1, block.KMAX-1)
            print('\t%d x %d x %d (%d)' %(block.IMAX, block.JMAX, block.KMAX, (block.IMAX-1)*(block.JMAX-1)*(nk)))

    status('Opening output exodusii file "%s"...' % output_file)
    exo = hpg.exodusii.exodusii_file(output_file, 'w')
    ndim = 3
    nnodes = 8
    type = 'HEX'
    if blocks[0].KMAX == 1:
        ndim = 2
        nnodes = 4
        type = 'QUAD'
    name = 'converted'
    node_count = blocks[0].IMAX * blocks[0].JMAX * blocks[0].KMAX
    nk = max(1, blocks[0].KMAX-1)
    cell_count = (blocks[0].IMAX - 1) * (blocks[0].JMAX - 1) * (nk)
    sidesets = [] # Maybe do something here
    exo.put_init(name, ndim, node_count, cell_count, len(blocks), 0, len(sidesets))
    status('Initialized exodusii file "%s".' % output_file)
    exo.put_coord(np.reshape(blocks[0].X, (node_count), order='F'),
                  np.reshape(blocks[0].Y, (node_count), order='F'),
                  np.reshape(blocks[0].Z, (node_count), order='F'))
    status('Wrote coordinates to exodusii file "%s".' % output_file)
    # Write out the blocks
    status('Wrote block(s) to exodusii file "%s".' % output_file)
    for nb, block in enumerate(blocks):
        id = nb + 1
        exo.put_element_block(id, type, cell_count, nnodes)
        exo.put_element_block_name(id, 'block_%d'%id)
        #
        conn = []
        if ndim == 2:
            k = 0
            for j in range(block.JMAX - 1):
                for i in range(block.IMAX - 1):
                    cell = [node_index(i,   j,   k,   block.IMAX, block.JMAX, block.KMAX),
                            node_index(i+1, j,   k,   block.IMAX, block.JMAX, block.KMAX),
                            node_index(i+1, j+1, k,   block.IMAX, block.JMAX, block.KMAX),
                            node_index(i,   j+1, k,   block.IMAX, block.JMAX, block.KMAX)]
                    conn.append(cell)
        else:
            for k in range(block.KMAX - 1):
                for j in range(block.JMAX - 1):
                    for i in range(block.IMAX - 1):
                        cell = [node_index(i,   j,   k,   block.IMAX, block.JMAX, block.KMAX),
                                node_index(i+1, j,   k,   block.IMAX, block.JMAX, block.KMAX),
                                node_index(i+1, j+1, k,   block.IMAX, block.JMAX, block.KMAX),
                                node_index(i,   j+1, k,   block.IMAX, block.JMAX, block.KMAX),
                                node_index(i,   j,   k+1, block.IMAX, block.JMAX, block.KMAX),
                                node_index(i+1, j,   k+1, block.IMAX, block.JMAX, block.KMAX),
                                node_index(i+1, j+1, k+1, block.IMAX, block.JMAX, block.KMAX),
                                node_index(i,   j+1, k+1, block.IMAX, block.JMAX, block.KMAX)]
                        conn.append(cell)
        exo.put_element_conn(id, np.array(conn))
    status('Done.')

    exo.close()
