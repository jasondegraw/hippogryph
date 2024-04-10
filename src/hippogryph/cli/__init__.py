# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import click
import os
import hippogryph as hpg
import enum

from ..__about__ import __version__

class Format(enum.Enum):
    EXO = enum.auto()
    PLOT3D = enum.auto()
    #PLOT3D_3DMB = enum.auto()
    #PLOT3D_3DSB = enum.auto()
    PLOT3D_2D = enum.auto()
    #PLOT3D_2DMB = enum.auto()
    #PLOT3D_2DSB = enum.auto()

def guess_output_format(filename: str) -> Format|None:
    name, extension = os.path.splitext(filename)
    match extension:
        case 'exo':
            return Format.EXO
        case 'xyz':
            return Format.PLOT3D
        case 'xy':
            return Format.PLOT3D_2D
    return None

@click.command()
@click.option('-x', '--x-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=32.0, help='Length of the grid in the x direction.')
@click.option('-y', '--y-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=1.0, help='Length of the grid in the y direction.')
#@click.option('-z', '--z-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the z direction.')
@click.option('-o', '--output', type=click.File('w'), show_default=True, default='hgf.exo', help='Write output to the specified file.')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d', 'plot3d_2d']), default=None, help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
def channel(x_length, y_length, output, format, ascii):
    '''
    Generate a channel grid.
    '''
    block = hpg.Block('domain')
    box = hpg.Box(ni=32, nj=32, block=block, left_label='west', right_label='east',
                  up_label='north', down_label='south')
    mesh = hpg.Mesh('channel')
    mesh.add(box)
    mesh.build()

    xgrid = hpg.Uniform.from_intervals(1.0, x_length)
    ygrid = hpg.Uniform.from_intervals(1.0, y_length, shift=-0.5)
    #xunif = hpg.Uniform.from_delta(ygrid.delta, 16*N)
    #xstretch = hpg.Geometric.from_delta(xunif.delta, 16, N)
    #xgrid = hpg.Composite([xunif, xstretch])

    mesh.apply(xgrid=xgrid, ygrid=ygrid)

def validate_even_int(ctx: click.core.Context,
                  param: click.core.Argument, value: str) -> int:
    try:
        v = int(value)
        if v <= 0:
            raise click.BadParameter('Number of elements must be positive.')
        if v % 2 != 0:
            raise click.BadParameter('Number of elements must be even.')
        return v
    except ValueError:
        raise click.BadParameter('Number of elements must be an integer.')

@click.command()
@click.option('-n', '--number', callback=validate_even_int, show_default=True, default=32, help='Number of elements across the channel (must be even).')
#@click.option('-z', '--z-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the z direction.')
@click.option('-o', '--output', type=click.Path(writable=True, dir_okay=False), show_default=True, default='hpg.exo', help='Write output to the specified file.')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d', 'plot3d_2d']), default=None, help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
def bfs(number, output, format, ascii):
    '''
    Generate a backward-facing step grid
    '''
    mesh = hpg.backward_step(int(number/2.0))
    mesh.write_exodusii(output)


@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='hippogryph')
@click.pass_context
def hippogryph(ctx: click.Context):
    pass

hippogryph.add_command(channel)
hippogryph.add_command(bfs)