# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
from . import exodusii
from . import plot3d
from .__about__ import __version__
import struct
import numpy as np
import textwrap

def discard_output(mesg):
    pass

def node_index(i, j, k, ni, nj, nk):
    return i + 1 + j * ni + k * ni * nj

def read_plotXd(input_file:str, status=discard_output):
    # Just assume it exists
    status('Attempting to read "%s" as a 3D, multiblock, unformatted file...' % input_file)
    try:
        return plot3d.read_plot3D(input_file)
    except struct.error:
        # Nope, not 3D, multiblock, unformatted
        pass
    status('Attempting to read "%s" as a 3D, multiblock, formatted file...' % input_file)
    try:
        return plot3d.read_plot3D(input_file, binary=False)
    except IndexError:
        # Not 3D, multiblock, formatted
        pass
    status('Attempting to read "%s" as a 2D, multiblock, unformatted file...' % input_file)
    try:
        return plot3d.read_plot2D(input_file)
    except struct.error:
        # Not 2D, multiblock, unformatted
        pass
    return plot3d.read_plot2D(input_file, binary=False)

def convert_plot3d(input_file, output_file, verbose=False):

    status = discard_output
    if verbose:
        status = print

    blocks = read_plotXd(input_file, status=status)

    if len(blocks) > 1:
        raise NotImplementedError('Conversion of multiblock grid is not yet supported')

    status('Read %d block(s)' % len(blocks))
    if verbose:
        for block in blocks:
            nk = max(1, block.KMAX-1)
            print('\t%d x %d x %d (%d)' %(block.IMAX, block.JMAX, block.KMAX, (block.IMAX-1)*(block.JMAX-1)*(nk)))

    status('Opening output exodusii file "%s"...' % output_file)
    exo = exodusii.exodusii_file(output_file, 'w')
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

def exo_info_lines(args:list[str], max_line_length:int=80, indent:int=2):
    offset = indent*' '
    result = [f'hippogryph v{__version__}']
    txt = ' '.join(args)
    return [f'hippogryph v{__version__}'] + textwrap.wrap(txt, width=max_line_length-1, break_on_hyphens=False,
                                                          subsequent_indent=offset)
    current = args[0]
    for arg in args[1:]:
        length = len(arg)
        if length >= max_line_length:
            pass
            #return []
        if len(current) + 1 + length >= max_line_length:
            result.append(current)
            current = offset + arg
        else:
            current += ' ' + arg
    return result
