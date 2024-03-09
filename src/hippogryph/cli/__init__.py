# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import click
import os
import hippogryph as hgf
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
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d', 'plot3d_2d']), default=True, default='exo', help='Write summary in JSON format.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
@click.option('-x', '--x-length', type=click.File('w'), show_default=True, default=32.0, help='Length of the grid in the x direction.')
@click.option('-y', '--y-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the y direction.')
@click.option('-z', '--z-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the z direction.')
@click.option('-o', '--output', type=click.File('w'), show_default=True, default='hgf.exo', help='Write output to the specified file.')
def channel(epjson, json, output):
    click.echo('Generate channel grid')

@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='hippogryph')
@click.pass_context
def hippogryph(ctx: click.Context):
    pass

hippogryph.add_command(channel)