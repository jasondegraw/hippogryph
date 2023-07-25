# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import numpy
import exodusii

class SubBlock:
    def __init__(self, number, i, j, k, ni, nj, nk):
        self.number = number
        self.i = i
        self.j = j
        self.k = k
        self.ni = ni
        self.nj = nj
        self.nk = nk


class Block:
    def __init__(self, name: str, id: int = 1):
        self.name = name
        self.id = id
        self.element_count = 0

class Box:
    def __init__(self, i=0, j=0, k=0, ni=1, nj=1, nk=None, name=None, block=None,
                 front_label=None, back_label=None, left_label=None, right_label=None,
                 up_label=None, down_label=None):
        self.name = name
        self.i = i
        self.j = j
        self.k = k
        self.ni = ni
        self.nj = nj
        self.nk = nk
        if nk is None:
            self.k = 0
            self.nk = 1
            self.two_dimensional = True
        elif k is None:
            self.k = 0
            self.nk = 1
            self.two_dimensional = True
        #self.two_dimensional = k is None or nk is None
        self.block = block
        self.subsets = {}
        self.front = front_label
        if self.front:
            self.add_to_subsets(self.front, 
                                SubBlock(6, self.i, self.j, self.k+self.nk, self.ni, self.nj, 1))
        self.back = back_label
        if self.back:
            self.add_to_subsets(self.back, 
                                SubBlock(5, self.i, self.j, self.k, self.ni, self.nj, 1))
        self.left = left_label
        if self.left:
            self.add_to_subsets(self.left,
                                SubBlock(4, self.i, self.j, self.k, 1, self.nj, self.nk))
        self.right = right_label
        if self.right:
            self.add_to_subsets(self.right, 
                                SubBlock(2, self.i+self.ni, self.j, self.k, 1, self.nj, self.nk))
        self.up = up_label
        if self.up:
            self.add_to_subsets(self.up,
                                SubBlock(3, self.i, self.j+self.nj, self.k, self.ni, 1, self.nk))
        self.down = down_label
        if self.down:
            self.add_to_subsets(self.down,
                                SubBlock(1, self.i, self.j, self.k, self.ni, 1, self.nk))   

    def add_to_subsets(self, name, obj):
        if name in self.subsets:
            self.subsets[name].append(obj)
        else:
            self.subsets[name] = [obj]

    def footprint(self):
        return numpy.ones((self.ni, self.nj, self.nk), dtype=numpy.uint8)
    
    def sidesets(self):
        return list(self.subsets.keys())

class AlreadyMeshed(Exception):
    pass

class Mesh:
    def __init__(self, name):
        self.name = name
        self._primitives = []
        self.two_dimensional = False
        self._meshed = False

    @classmethod
    def from_array(cls, name, shape, array): # rework this with Numpy or something
        primitives = []
        if len(shape) == 1:
            raise NotImplementedError
        elif len(shape) == 2:
            index = 0
            imax = 0
            jmax = 0
            for j in range(shape[1]):
                last_j = jmax
                row = []
                for i in range(shape[0]):
                    row.append(array[index])
                    jmax = max(jmax, array[index].nj)
                    index += 1
                for obj in row:
                    obj.i += imax
                    obj.j += last_j
                    imax = obj.i
                    if obj.nj != jmax:
                        raise NotImplementedError
                last_j += jmax
                primitives.extend(row)
        elif len(shape) > 2:
            raise NotImplementedError
        object = cls(name)
        for primitive in primitives:
            object.add(primitive)

        return object

    @property
    def primitives(self):
        return self._primitives
    
    def add(self, primitive):
        if self._primitives:
            if primitive.two_dimensional == self.two_dimensional:
                self._primitives.append(primitive)
            else:
                return False
        else:
            self._primitives.append(primitive)
            self.two_dimensional = primitive.two_dimensional

    def mesh(self, force=False):
        if self._meshed:
            raise AlreadyMeshed('Mesh "{self.name}" is already meshed')
        if force:
            raise NotImplementedError
        
        self.build(force=force)
    
    def apply(self, xgrid, ygrid, zgrid=None, force=False):
        if self._meshed:
            raise AlreadyMeshed('Mesh "{self.name}" is already meshed')
        if force:
            raise NotImplementedError
        
        if self.two_dimensional:
            x = numpy.zeros(self.ni+1)
            for i in range(self.ni+1):
                x[i] = xgrid.s(i)
            y = numpy.zeros(self.nj+1)
            for j in range(self.nj+1):
                y[j] = ygrid.s(j)

            self.x = numpy.zeros(self.node_count)
            self.y = numpy.zeros(self.node_count)

            k = 0
            index = 0
            for j in range(self.nj+1):
                for i in range(self.ni+1):
                    if self.node_index[i, j, k] > 0:
                        self.x[index] = x[i]
                        self.y[index] = y[j]
                        index += 1
        else:
            raise NotImplementedError

    def build(self, force=False):
        if self._meshed:
            raise AlreadyMeshed('Mesh "{self.name}" is already meshed')
        if force:
            raise NotImplementedError
        
        if not self._primitives:
            return
        
        # Get the blocks and sidesets we have
        self.blocks = []
        sidesets = set()
        for primitive in self._primitives:
            if not primitive.block in self.blocks:
                self.blocks.append(primitive.block)
            sidesets.update(primitive.sidesets())
            print(primitive.sidesets())
        
        # Number the blocks from 1
        reverse_lookup = {}
        for i, block in enumerate(self.blocks):
            block.id = i+1
            reverse_lookup[block.id] = block

        # Get the sidesets we have
        
        # Figure out the size
        if self.two_dimensional:
            i_offset = self._primitives[0].i
            j_offset = self._primitives[0].j
            i_max = self._primitives[0].i + self._primitives[0].ni
            j_max = self._primitives[0].j + self._primitives[0].nj
            for block in self._primitives[1:]:
                i_offset = min(i_offset, block.i)
                j_offset = min(j_offset, block.j)
                i_max = max(i_max, block.i + block.ni)
                j_max = max(j_max, block.j + block.nj)
            self.ni = i_max - i_offset
            self.nj = j_max - j_offset
            self.nk = 1
            self.i_offset = i_offset
            self.j_offset = j_offset
            self.k_offset = 0
        else:
            i_offset = self._primitives[0].i
            j_offset = self._primitives[0].j
            k_offset = self._primitives[0].k
            i_max = self._primitives[0].i + self._primitives[0].ni
            j_max = self._primitives[0].j + self._primitives[0].nj
            k_max = self._primitives[0].k + self._primitives[0].nk
            for block in self._primitives[1:]:
                i_offset = min(i_offset, block.i)
                j_offset = min(j_offset, block.j)
                k_offset = min(k_offset, block.k)
                i_max = max(i_max, block.i + block.ni)
                j_max = max(j_max, block.j + block.nj)
                k_max = max(k_max, block.k + block.nk)
            self.ni = i_max - i_offset
            self.nj = j_max - j_offset
            self.nk = k_max - k_offset
            self.i_offset = i_offset
            self.j_offset = j_offset
            self.k_offset = k_offset

        self.cells = numpy.zeros((self.ni, self.nj, self.nk), dtype=numpy.uint8)
        self.cell_index = numpy.zeros((self.ni, self.nj, self.nk), dtype=numpy.uint64)
        self.node_index = numpy.zeros((self.ni+1, self.nj+1, self.nk+1), dtype=numpy.uint64)

        # Map it out
        for primitive in self.primitives:
            k0 = 0
            k1 = 1
            if not self.two_dimensional:
                k0 = primitive.k - self.k_offset
                k1 = primitive.k + primitive.nk - self.k_offset
            self.cells[primitive.i - self.i_offset:primitive.i + primitive.ni - self.i_offset,
                       primitive.j - self.j_offset:primitive.j + primitive.nj - self.j_offset,
                       k0:k1] = primitive.block.id * primitive.footprint()
        
        # Assign cell numbers and flag nodes, probably need to do this differently
        index = 0
        if self.two_dimensional:
            k = 0
            for j in range(self.nj):
                for i in range(self.ni):
                    if self.cells[i, j, k] > 0:
                        index += 1
                        self.cell_index[i, j, k] = index
                        self.node_index[i,   j,   k] = 1
                        self.node_index[i+1, j,   k] = 1
                        self.node_index[i+1, j+1, k] = 1
                        self.node_index[i,   j+1, k] = 1
                        reverse_lookup[self.cells[i, j, k]].element_count += 1
        else:
            for k in range(self.nk):
                for j in range(self.nj):
                    for i in range(self.ni):
                        if self.cells[i, j, k] > 0:
                            index += 1
                            self.cell_index[i, j, k] = index
                            self.node_index[i,   j,   k] = 1
                            self.node_index[i+1, j,   k] = 1
                            self.node_index[i+1, j+1, k] = 1
                            self.node_index[i,   j+1, k] = 1
                            self.node_index[i,   j,   k+1] = 1
                            self.node_index[i+1, j,   k+1] = 1
                            self.node_index[i+1, j+1, k+1] = 1
                            self.node_index[i,   j+1, k+1] = 1
                            reverse_lookup[self.cells[i, j, k]].element_count += 1
        self.cell_count = index
        # Now number the nodes
        k0 = self.nk + 1
        if self.two_dimensional:
            k0 = 1
        index = 0
        for k in range(k0):
            for j in range(self.nj+1):
                for i in range(self.ni+1):
                    if self.node_index[i, j, k] > 0:
                        index += 1
                        self.node_index[i, j, k] = index
        self.node_count = index
        print(self.cell_count, self.node_count)
        print(len(sidesets))
        print(sidesets)
        print(self.ni, self.nj)

    def save(self, filename):
        exo = exodusii.exodusii_file(filename, 'w')
        ndim = 3
        nnodes = 8
        type = 'HEX'
        if self.two_dimensional:
            ndim = 2
            nnodes = 4
            type = 'QUAD'
        sidesets = []
        exo.put_init(self.name, ndim, self.node_count, self.cell_count,
                     len(self.blocks), 0, len(sidesets))
        exo.put_coord(self.x, self.y)

        for nb, block in enumerate(self.blocks):
            id = nb + 1
            exo.put_element_block(id, type, block.element_count, nnodes)
            exo.put_element_block_name(id, block.name)
            #
            conn = []
            if self.two_dimensional:
                k = 0
                for j in range(self.nj):
                    for i in range(self.ni):
                        if self.cells[i, j, k] == block.id:
                            cell = [self.node_index[i,   j,   k],
                                    self.node_index[i+1, j,   k],
                                    self.node_index[i+1, j+1, k],
                                    self.node_index[i,   j+1, k]]
                            conn.append(cell)
            else:
                for k in range(self.nk):
                    for j in range(self.nj):
                        for i in range(self.ni):
                            if self.cells[i, j, k] > 0:
                                cell = [self.node_index[i,   j,   k],
                                        self.node_index[i+1, j,   k],
                                        self.node_index[i+1, j+1, k],
                                        self.node_index[i,   j+1, k],
                                        self.node_index[i,   j,   k+1],
                                        self.node_index[i+1, j,   k+1],
                                        self.node_index[i+1, j+1, k+1],
                                        self.node_index[i,   j+1, k+1]]
                                conn.append(cell)
            exo.put_element_conn(block.id, numpy.array(conn))

        exo.close()
