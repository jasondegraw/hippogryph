# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import hippogryph
import numpy as np
import os
import hippogryph.plot3d

def dev_null(mesg):
    pass

pts = [[0.0, -0.75], [0.5, -0.75], [1.0, -0.75], [1.5, -0.75], [2.0, -0.75],
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

def test_bfs(tmpdir):
    mesh = hippogryph.backward_step(1)
    assert len(mesh.blocks) == 1
    ib = mesh.blocks[0].iblank()
    assert ib.shape == (35,3,1)
    np.testing.assert_array_equal(ib, 1)
    #fullpath = os.path.join(tmpdir,'bfs.exo')
    #mesh.write_exodusii(fullpath)
    #with exodusii.File(fullpath, mode="r") as exof:
    #    domain = exof.get_element_block(1)
    #    assert domain.name == 'domain'
    #    coords = exof.get_coords()
    #    assert np.allclose(coords, pts)
    #    assert len(coords) == 105
    #    conn = exof.get_element_connectivity(1).elem_conn
    #    assert len(conn) == 68

tee_ib_2 = [[1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]]

tee_ib_3 = np.expand_dims(tee_ib_2, axis=2)

def test_tee(tmpdir):
    mesh = hippogryph.tee_junction(1)
    assert len(mesh.blocks) == 1
    ib = mesh.blocks[0].iblank()
    assert ib.shape == (13,10,1)
    assert tee_ib_3.shape == (13,10,1)
    np.testing.assert_array_equal(ib, tee_ib_3)
    #fullpath = os.path.join(tmpdir,'bfs.exo')
    #mesh.write_exodusii(fullpath)
    #with exodusii.File(fullpath, mode="r") as exof:
    #    domain = exof.get_element_block(1)
    #    assert domain.name == 'domain'
    #    coords = exof.get_coords()
    #    assert np.allclose(coords, pts)
    #    assert len(coords) == 105
    #    conn = exof.get_element_connectivity(1).elem_conn
    #    assert len(conn) == 68

x = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 
              0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0,
              0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0])

y = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
              -0.5, -0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
              -0.5, -0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5])

z = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
              0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

def test_channel_3d_primitive(tmpdir):
    mesh = hippogryph.channel(x=2.0, z=1.0, ni=4, nj=2, nk=2)
    assert len(mesh.blocks) == 1
    ib = mesh.blocks[0].iblank()
    assert ib.shape == (5, 3, 3)
    assert np.allclose(ib, 1)
    assert mesh.node_count == 45
    assert mesh.cell_count == 16
    assert mesh.blocks[0].x.shape == (45,)
    assert np.allclose(mesh.blocks[0].x, x)
    assert mesh.blocks[0].y.shape == (45,)
    assert np.allclose(mesh.blocks[0].y, y)
    assert mesh.blocks[0].z.shape == (45,)
    assert np.allclose(mesh.blocks[0].z, z)
    block = hippogryph.plot3d.Block(np.reshape(mesh.blocks[0].x, (5,3,3), order='F'),
                                    np.reshape(mesh.blocks[0].y, (5,3,3), order='F'),
                                    np.reshape(mesh.blocks[0].z, (5,3,3), order='F'))
    filename = os.path.join(tmpdir, 'out.xyz')
    hippogryph.plot3d.write_plot3D(filename, [block])
    blocks = hippogryph.plot3d.read_plot3D(filename) #, binary:bool=True,big_endian:bool=False,read_double:bool=True):
    assert len(blocks) == 1
    assert blocks[0].X.shape == (5,3,3)
    assert blocks[0].Y.shape == (5,3,3)
    assert blocks[0].Z.shape == (5,3,3)
    assert np.allclose(np.reshape(blocks[0].X, (45,), order='F'), x)
    assert np.allclose(np.reshape(blocks[0].Y, (45,), order='F'), y)
    assert np.allclose(np.reshape(blocks[0].Z, (45,), order='F'), z)

def test_channel_3d(tmpdir):
    mesh = hippogryph.channel(x=2.0, z=1.0, ni=4, nj=2, nk=2)
    assert len(mesh.blocks) == 1
    ib = mesh.blocks[0].iblank()
    assert ib.shape == (5, 3, 3)
    assert np.allclose(ib, 1)
    assert mesh.node_count == 45
    assert mesh.cell_count == 16
    assert mesh.blocks[0].x.shape == (45,)
    assert np.allclose(mesh.blocks[0].x, x)
    assert mesh.blocks[0].y.shape == (45,)
    assert np.allclose(mesh.blocks[0].y, y)
    assert mesh.blocks[0].z.shape == (45,)
    assert np.allclose(mesh.blocks[0].z, z)
    filename = os.path.join(tmpdir, 'chan.xyz')
    mesh.write_plot3d(filename)
    blocks = hippogryph.plot3d.read_plot3D(filename) #, binary:bool=True,big_endian:bool=False,read_double:bool=True):
    assert len(blocks) == 1
    assert blocks[0].X.shape == (5,3,3)
    assert blocks[0].Y.shape == (5,3,3)
    assert blocks[0].Z.shape == (5,3,3)
    assert np.allclose(np.reshape(blocks[0].X, (45,), order='F'), x)
    assert np.allclose(np.reshape(blocks[0].Y, (45,), order='F'), y)
    assert np.allclose(np.reshape(blocks[0].Z, (45,), order='F'), z)