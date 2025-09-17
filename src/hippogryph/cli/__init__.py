# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import click
import sys
import hippogryph as hpg

from ..__about__ import __version__

@click.command()
@click.option('-x', '--x-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=32.0, help='Length of the grid in the x direction.')
@click.option('-y', '--y-length', type=click.FloatRange(0.0, min_open=True), show_default=True, default=1.0, help='Length of the grid in the y direction.')
@click.option('-z', '--z-length', type=click.FloatRange(0.0, min_open=False), show_default=True, default=0.0, help='Length of the grid in the z direction.')
@click.option('--dy', type=click.FloatRange(0.0, min_open=True), show_default=True, default=None, help='Length of first grid spacing in the j (y) direction.')
@click.option('--yplus', type=click.FloatRange(0.0, min_open=True), show_default=True, default=None, help='y+ of the first grid spacing in the j (y) direction.')
@click.option('--re-tau', type=click.FloatRange(0.0, min_open=True), show_default=True, default=None, help='Friction Reynolds number to use for the j (y) direction spacings.')
@click.option('-i', '--ni', type=click.IntRange(1), show_default=True, default=32, help='Number of cells in the i (x) direction.')
@click.option('-j', '--nj', type=click.IntRange(1), show_default=True, default=32, help='Number of cells in the j (y) direction.')
@click.option('-k', '--nk', type=click.IntRange(0), show_default=True, default=0, help='Number of cells in the k (z) direction.')
@click.option('-o', '--output', type=click.Path(dir_okay=False, writable=True), show_default=True, default=None, help='File to write output to, defaults to "halfchan.exo|xyz|xy".')
@click.option('-f', '--format', type=click.Choice(['exo', 'plot3d']), default='exo', help='Specify format to use.')
@click.option('-a', '--ascii', is_flag=True, show_default=True, default=False, help='Write ASCII format (if possible).')
@click.option('-w', '--wall-label', type=str, show_default=True, default='wall', help='Label for the wall boundary.')
@click.option('-c', '--centerline-label', type=str, show_default=True, default='centerline', help='Label for the centerline boundary.')
@click.option('-l', '--left-label', type=str, show_default=True, default='inflow', help='Label for the left boundary.')
@click.option('-r', '--right-label', type=str, show_default=True, default='outflow', help='Label for the right boundary.')
@click.option('-n', '--name', type=str, show_default=True, default='domain', help='Label for the flow domain.')
def half_channel(x_length, y_length, z_length, dy, yplus, re_tau, ni, nj, nk, output, format, ascii, wall_label,
                 centerline_label, left_label, right_label, name):
    """
    Generate a boundary layer grid.
    """
    binary = not ascii
    if yplus is not None and re_tau is not None:
        # Put y+ = Re_\tau at the y=y_length boundary:
        #   u_tau * y_length/ nu = re_tau => u_tau / nu = re_tau / y_length
        # and then
        #   yplus = dy * u_tau / nu => dy = yplus * y_length / re_tau
        dy = yplus * y_length / re_tau
    mesh = hpg.half_channel(x=x_length, y=y_length, z=z_length, ni=ni, nj=nj, nk=nk, dy=dy, left_label=left_label,
                            right_label=right_label, down_label=wall_label, up_label=centerline_label, domain_label=name)
    if format == 'exo':
        if output is None:
            output = 'halfchan.exo'
        success = mesh.write_exodusii(output, info=hpg.exo_info_lines(sys.argv))
    elif format == 'plot3d':
        if output is None:
            output = 'halfchan.xyz'
            if mesh.two_dimensional:
                output = 'halfchan.xy'
        success = mesh.write_plot3d(output, binary=binary)
    if not success:
        click.echo('Writing output to "%s" failed' % output)
    click.echo('# Mesh Statistics #')
    click.echo('x extents: %e to %e' % (mesh.blocks[0].x[0], mesh.blocks[0].x[-1]))
    click.echo('y extents: %e to %e' % (mesh.blocks[0].y[0], mesh.blocks[0].y[-1]))

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
@click.option('-t', '--top-wall-label', type=str, show_default=True, default='wall', help='Label for the wall boundary.')
@click.option('-b', '--bottom-wall-label', type=str, show_default=True, default='centerline', help='Label for the centerline boundary.')
@click.option('-l', '--left-label', type=str, show_default=True, default='inflow', help='Label for the left boundary.')
@click.option('-r', '--right-label', type=str, show_default=True, default='outflow', help='Label for the right boundary.')
@click.option('-n', '--name', type=str, show_default=True, default='domain', help='Label for the flow domain.')
def channel(x_length, y_length, z_length, ni, nj, nk, output, format, ascii, top_wall_label, bottom_wall_label,
            left_label, right_label, name):
    """
    Generate a channel grid.
    """
    binary = not ascii
    mesh = hpg.channel(x=x_length, y=y_length, z=z_length, ni=ni, nj=nj, nk=nk, left_label=left_label, right_label=right_label,
                       up_label=top_wall_label, down_label=bottom_wall_label, domain_label=name)
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
        click.echo('Writing output to "%s" failed' % output)

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
        click.echo('Writing output to "%s" failed' % output)

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
        click.echo('Writing output to "%s" failed' % output)

@click.command()
@click.option('-o', '--output', type=click.Path(writable=True, dir_okay=False), show_default=True, default='p3d.exo', help='File to write output to, defaults to "p3d.exo".')
@click.option('-v', '--verbose', is_flag=True, show_default=True, default=False, help='Operate verbosely and write out progress information.')
@click.argument('plot3d-file', type=click.Path(exists=True))
def convert_plot3d(output, verbose, plot3d_file):
    '''
    Convert a Plot3D grid to Exodus II
    '''
    hpg.convert_plot3d(plot3d_file, output, verbose=verbose)


@click.group(context_settings={'help_option_names': ['-h', '--help']}, invoke_without_command=False)
@click.version_option(version=__version__, prog_name='hippogryph')
@click.pass_context
def hippogryph(ctx: click.Context):
    pass

hippogryph.add_command(half_channel)
hippogryph.add_command(channel)
hippogryph.add_command(bfs)
hippogryph.add_command(tjunct)
hippogryph.add_command(convert_plot3d)