# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import hippogryph
import hippogryph.exodusii as exodusii
import numpy as np
import os

import hippogryph.plot3d

def dev_null(mesg):
    pass

bfs_pts = [[0.0, -0.75], [0.5, -0.75], [1.0, -0.75], [1.5, -0.75], [2.0, -0.75],
           [2.5, -0.75], [3.0, -0.75], [3.5, -0.75], [4.0, -0.75], [4.5, -0.75],
           [5.0, -0.75], [5.5, -0.75], [6.0, -0.75], [6.5, -0.75], [7.0, -0.75],
           [7.5, -0.75], [8.0, -0.75], [8.5, -0.75], [9.0, -0.75], [9.5, -0.75],
           [10.0, -0.75], [10.5, -0.75], [11.0, -0.75], [11.5, -0.75], [12.0, -0.75],
           [12.5, -0.75], [13.0, -0.75], [13.5, -0.75], [14.0, -0.75], [14.5, -0.75],
           [15.0, -0.75], [15.5, -0.75], [16.0, -0.75], [16.5, -0.75], [32.0, -0.75],
           [0.0, -0.25], [0.5, -0.25], [1.0, -0.25], [1.5, -0.25], [2.0, -0.25],
           [2.5, -0.25], [3.0, -0.25], [3.5, -0.25], [4.0, -0.25], [4.5, -0.25],
           [5.0, -0.25], [5.5, -0.25], [6.0, -0.25], [6.5, -0.25], [7.0, -0.25],
           [7.5, -0.25], [8.0, -0.25], [8.5, -0.25], [9.0, -0.25], [9.5, -0.25],
           [10.0, -0.25], [10.5, -0.25], [11.0, -0.25], [11.5, -0.25], [12.0, -0.25],
           [12.5, -0.25], [13.0, -0.25], [13.5, -0.25], [14.0, -0.25], [14.5, -0.25],
           [15.0, -0.25], [15.5, -0.25], [16.0, -0.25], [16.5, -0.25], [32.0, -0.25],
           [0.0, 0.25], [0.5, 0.25], [1.0, 0.25], [1.5, 0.25], [2.0, 0.25],
           [2.5, 0.25], [3.0, 0.25], [3.5, 0.25], [4.0, 0.25], [4.5, 0.25],
           [5.0, 0.25], [5.5, 0.25], [6.0, 0.25], [6.5, 0.25], [7.0, 0.25],
           [7.5, 0.25], [8.0, 0.25], [8.5, 0.25], [9.0, 0.25], [9.5, 0.25],
           [10.0, 0.25], [10.5, 0.25], [11.0, 0.25], [11.5, 0.25], [12.0, 0.25],
           [12.5, 0.25], [13.0, 0.25], [13.5, 0.25], [14.0, 0.25], [14.5, 0.25],
           [15.0, 0.25], [15.5, 0.25], [16.0, 0.25], [16.5, 0.25], [32.0, 0.25]]

tee_pts = [[0.0, -0.5], [1.0, -0.5], [2.0, -0.5], [3.0, -0.5], [4.0, -0.5],
           [5.0, -0.5], [6.0, -0.5], [7.0, -0.5], [8.0, -0.5], [9.0, -0.5],
           [10.0, -0.5], [11.0, -0.5], [18.0, -0.5],
           [0.0, 0.5], [1.0, 0.5], [2.0, 0.5], [3.0, 0.5], [4.0, 0.5],
           [5.0, 0.5], [6.0, 0.5], [7.0, 0.5], [8.0, 0.5], [9.0, 0.5],
           [10.0, 0.5], [11.0, 0.5], [18.0, 0.5],
           [3.0, 1.5], [4.0, 1.5],
           [3.0, 2.5], [4.0, 2.5],
           [3.0, 3.5], [4.0, 3.5],
           [3.0, 4.5], [4.0, 4.5],
           [3.0, 5.5], [4.0, 5.5],
           [3.0, 6.5], [4.0, 6.5],
           [3.0, 7.5], [4.0, 7.5],
           [3.0, 14.5], [4.0, 14.5]]

channel_pts = [[0.0, -0.5], [0.5, -0.5], [1.0, -0.5], [1.5, -0.5], [2.0, -0.5], [2.5, -0.5], [3.0, -0.5], [3.5, -0.5], [4.0, -0.5],
               [0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.5, 0.0], [2.0, 0.0], [2.5, 0.0], [3.0, 0.0], [3.5, 0.0], [4.0, 0.0],
               [0.0, 0.5], [0.5, 0.5], [1.0, 0.5], [1.5, 0.5], [2.0, 0.5], [2.5, 0.5], [3.0, 0.5], [3.5, 0.5], [4.0, 0.5],]

def test_bfs_2d(tmpdir):
    mesh = hippogryph.backward_step(1)
    assert len(mesh.blocks) == 1
    fullpath = os.path.join(tmpdir,'bfs2d.exo')
    mesh.write_exodusii(fullpath)
    with exodusii.File(fullpath, mode="r") as exof:
        domain = exof.get_element_block(1)
        assert domain.name == 'domain'
        coords = exof.get_coords()
        assert len(coords) == 105
        assert np.allclose(coords, bfs_pts)
        conn = exof.get_element_connectivity(1).elem_conn
        assert len(conn) == 68
        
def test_tee_2d(tmpdir):
    mesh = hippogryph.tee_junction(1)
    fullpath = os.path.join(tmpdir,'tee2d.exo')
    mesh.write_exodusii(fullpath)
    with exodusii.File(fullpath, mode="r") as exof:
        domain = exof.get_element_block(1)
        assert domain.name == 'domain'
        coords = exof.get_coords()
        assert len(coords) == len(tee_pts)
        assert np.allclose(coords, tee_pts)
        conn = exof.get_element_connectivity(1).elem_conn
        assert len(conn) == 20

def test_channel_2d(tmpdir):
    mesh = hippogryph.channel(x=4.0, y=1.0, ni=8, nj=2, nk=0)
    assert mesh.node_count == 27
    assert mesh.cell_count == 16
    fullpath = os.path.join(tmpdir,'chan2d.exo')
    mesh.write_exodusii(fullpath)
    with exodusii.File(fullpath, mode="r") as exof:
        domain = exof.get_element_block(1)
        assert domain.name == 'domain'
        coords = exof.get_coords()
        assert len(coords) == 27
        assert np.allclose(coords, channel_pts)
        conn = exof.get_element_connectivity(1).elem_conn
        assert len(conn) == 16

def test_channel_3d(tmpdir):
    mesh = hippogryph.channel(x=1.0, y=1.0, z=1.0, ni=4, nj=8, nk=2)
    assert mesh.node_count == 5*9*3
    assert mesh.cell_count == 64
    fullpath = os.path.join(tmpdir,'chan3d.exo')
    mesh.write_exodusii(fullpath)
    