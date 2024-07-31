# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import click
import hippogryph as hpg

from ..__about__ import __version__

@click.command()
@click.option('-x', '--x-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=32.0, help='Length of the grid in the x direction.')
@click.option('-y', '--y-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=1.0, help='Length of the grid in the y direction.')
@click.option('-z', '--z-length', type=click.FloatRange(0.0, min_open=False), show_default=True, default=0.0, help='Length of the grid in the z direction.')
@click.option('-i', '--ni', type=click.IntRange(1), show_default=True, default=32, help='Number of cells in the i (x) direction.')
@click.option('-j', '--nj', type=click.IntRange(1), show_default=True, default=32, help='Number of cells in the j (y) direction.')
@click.option('-k', '--nk', type=click.IntRange(0), show_default=True, default=0, help='Number of cells in the k (z) direction.')
@click.option('-o', '--output', type=click.Path(dir_okay=False, writable=True), show_default=True, default=None, help='File to write output to, defaults to "chan.exo|xyz|xy".')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d']), default='exo', help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
def channel(x_length, y_length, z_length, ni, nj, nk, output, format, ascii):
    """
    Generate a channel grid.
    """
    binary = not ascii
    mesh = hpg.channel(x=x_length, y=y_length, z=z_length, ni=ni, nj=nj, nk=nk)
    if format == 'exo':
        if output is None:
            output = 'chan.exo'
        success = mesh.write_exodusii(output)
    elif format == 'plot3d':
        if output is None:
            output = 'chan.xyz'
            if mesh.two_dimensional:
                output = 'chan.xy'
        success = mesh.write_plot3d(output, binary=binary)
    if not success:
        print('Writing output to "%s" failed' % output)

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
@click.option('-o', '--output', type=click.Path(writable=True, dir_okay=False), show_default=True, default=None, help='File to write output to, defaults to "bfs.exo|xyz|xy".')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d']), default=None, help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
def bfs(number, output, format, ascii):
    '''
    Generate a backward-facing step grid
    '''
    binary = not ascii
    mesh = hpg.backward_step(int(number/2.0))
    if format == 'exo':
        if output is None:
            output = 'bfs.exo'
        success = mesh.write_exodusii(output)
    elif format == 'plot3d':
        if output is None:
            output = 'bfs.xyz'
            if mesh.two_dimensional:
                output = 'bfs.xy'
        success = mesh.write_plot3d(output, binary=binary)
    if not success:
        print('Writing output to "%s" failed' % output)

@click.command()
@click.option('-n', '--number', callback=validate_even_int, show_default=True, default=32, help='Number of elements across the channel (must be even).')
#@click.option('-z', '--z-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the z direction.')
@click.option('-o', '--output', type=click.Path(writable=True, dir_okay=False), show_default=True, default=None, help='File to write output to, defaults to "tjunct.exo|xyz|xy".')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d']), default=None, help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
def tjunct(number, output, format, ascii):
    '''
    Generate a tee-junction grid
    '''
    binary = not ascii
    mesh = hpg.tee_junction(number)
    if format == 'exo':
        if output is None:
            output = 'tjunct.exo'
        success = mesh.write_exodusii(output)
    elif format == 'plot3d':
        if output is None:
            output = 'tjunct.xyz'
            if mesh.two_dimensional:
                output = 'tjunct.xy'
        success = mesh.write_plot3d(output, binary=binary)
    if not success:
        print('Writing output to "%s" failed' % output)

@click.command()

#@click.option('-z', '--z-length', type=click.File('w'), show_default=True, default=1.0, help='Length of the grid in the z direction.')
@click.option('-o', '--output', type=click.Path(writable=True, dir_okay=False), show_default=True, default=None, help='File to write output to, defaults to "p3d.exo|xyz|xy".')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d']), default=None, help='Specify format to use.')
def convert_plot3d(number, output, format, ascii):
    '''
    Generate a tee-junction grid
    '''


@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='hippogryph')
@click.pass_context
def hippogryph(ctx: click.Context):
    pass

hippogryph.add_command(channel)
hippogryph.add_command(bfs)
hippogryph.add_command(tjunct)