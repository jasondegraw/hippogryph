# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import numpy as np
from . import exodusii
from . import plot3d

class SideSet:
    """An object representing the Exodus II side set.
    
    This object specifies a range of elements and the side number that makes up an Exodus II side set.

    Parameters
    ----------
    number:
        The Exodus II side number.
    i,j,k:
        The origin of the side set.
    ni,nj,nkL
        The extents of the side set. One of these should be 1.
    """
    def __init__(self, number:int, i:int, j:int, k:int, ni:int, nj:int, nk:int):
        self.number = number
        self.i = i
        self.j = j
        self.k = k
        self.ni = ni
        self.nj = nj
        self.nk = nk

class ElementSet:
    """A named grouping of elements.

    This object (not to be confused with Block) allows elements to be grouped together,
    corresponds to the Exodus II block.

    Parameters
    ----------
    name:
        The name of the element block.

    id: optional
        The Exodus II id of the block, defaults to 1.
    """
    def __init__(self, name: str, id: int = 1):
        self.name = name
        self.id = id
        self.element_count = 0

class Box:
    """A primitive that selects an box-shaped area for meshing.
    
    The most basic primitive, this object selects a box shaped area for meshing and associates it
    with an element set. 

    Parameters
    ----------
    i: optional
        The i-direction origin of the box, defaults to 0.
    j: optional
        The j-direction origin of the box, defaults to 0.
    k: optional
        The k-direction origin of the box, defaults to None for a two-dimensional box.
    ni: optional
        The number of i-direction intervals, defaults to 1.
    nj: optional
        The number of j direction intervals, defaults to 1.
    nk: optional
        The number of k-direction intervals, defaults to None for a two-dimensional box.
    name: optional
        The name of the box.
    element_set: optional
        The name of an element set to associate with, an element set object to be associated with, or None.
    back_label: optional
        The name of the side set to be created on the k boundary of the box. Defaults to None.
    left_label: optional
        The name of the side set to be created on the i boundary of the box. Defaults to None.
    right_label: optional
        The name of the side set to be created on the i+ni boundary of the box. Defaults to None.
    up_label: optional
        The name of the side set to be created on the j+nj boundary of the box. Defaults to None.
    down_label: optional
        The name of the side set to be created on the j boundary of the box. Defaults to None.

    Notes
    -----
    Giving more than one of the label parameters the same name will result in a combined sideset. For example,
    if down_label and right_label are the same, the eventual output sideset will include both of those boundaries.
    """
    def __init__(self, i:int=0, j:int=0, k:int=0, ni:int=1, nj:int=1, nk:int|None=None, name:str|None=None, 
                 element_set:ElementSet|str|None=None, front_label:str|None=None, back_label:str|None=None,
                 left_label:str|None=None, right_label:str|None=None, up_label:str|None=None,
                 down_label:str|None=None):
        self.name = name
        self.i = i
        self.j = j
        self.k = k
        self.ni = ni
        self.nj = nj
        self.nk = nk
        self.two_dimensional = False
        if nk is None or nk == 0:
            self.k = 0
            self.nk = 1
            self.two_dimensional = True
        elif k is None:
            self.k = 0
            self.nk = 1
            self.two_dimensional = True
        self.element_set = element_set
        #self.subsets = {}
        self.front = front_label
        #if self.front:
        #    self.add_to_subsets(self.front, 
        #                        SideSet(6, self.i, self.j, self.k+self.nk, self.ni, self.nj, 1))
        self.back = back_label
        #if self.back:
        #    self.add_to_subsets(self.back, 
        #                        SideSet(5, self.i, self.j, self.k, self.ni, self.nj, 1))
        self.left = left_label
        #if self.left:
        #    self.add_to_subsets(self.left,
        #                        SideSet(4, self.i, self.j, self.k, 1, self.nj, self.nk))
        self.right = right_label
        #if self.right:
        #    self.add_to_subsets(self.right, 
        #                        SideSet(2, self.i+self.ni, self.j, self.k, 1, self.nj, self.nk))
        self.up = up_label
        #if self.up:
        #    self.add_to_subsets(self.up,
        #                        SideSet(3, self.i, self.j + self.nj, self.k, self.ni, 1, self.nk))
        self.down = down_label
        #if self.down:
        #    self.add_to_subsets(self.down,
        #                        SideSet(1, self.i, self.j, self.k, self.ni, 1, self.nk))
            
        names = [self.down, self.right, self.up, self.left, self.back, self.front]
        numbers = [1, 2, 3, 4, 5, 6]
        self._sidesets = [(label, number) for label, number in zip(names, numbers) if label is not None]

    #def add_to_subsets(self, name, obj):
    #    if name in self.subsets:
    #        self.subsets[name].append(obj)
    #    else:
    #        self.subsets[name] = [obj]

    def footprint(self)->np.ndarray[np.uint8]:
        """Return a footprint array of 1s and 0s.
        
        This method returns an array of the correct shape for this primitive with each active node
        marked with a 1.
        
        Returns
        -------
        np.ndarray:
            The footprint array.
        """
        return np.ones((self.ni, self.nj, self.nk), dtype=np.uint8)
    
    def sidesets(self)->list[str]:
        """Returns the names of all of the sidesets.
        
        Returns
        -------
        list[str]
            The list of sideset names.
        """
        return [el[0] for el in self._sidesets]
    
    def sideset(self, name:str)->list[SideSet]:
        """Return the named sideset.
        
        This method will return a list of SideSet objects that are grouped under the given name.

        Parameters
        ----------
        name:
            The name of the sidesets to return.
        
        Returns
        -------
        list[SideSet]:
            The list of SideSet objects that have the given name.
        """
        numbers = [el[1] for el in self._sidesets if el[0] == name]
        results = []
        for number in numbers:
            if number == 1:
                results.append(SideSet(1, self.i, self.j, self.k, self.ni, 1, self.nk))
            elif number == 2:
                results.append(SideSet(2, self.i + self.ni - 1, self.j, self.k, 1, self.nj, self.nk))
            elif number == 3:
                results.append(SideSet(3, self.i, self.j + self.nj - 1, self.k, self.ni, 1, self.nk))
            elif number == 4:
                results.append(SideSet(4, self.i, self.j, self.k, 1, self.nj, self.nk))
            elif number == 5:
                results.append(SideSet(5, self.i, self.j, self.k, self.ni, self.nj, 1))
            else:
                results.append(SideSet(6, self.i, self.j, self.k + self.nk - 1, self.ni, self.nj, 1))
        return results

class AlreadyMeshed(Exception):
    """Raised if a mesh object has already been meshed."""
    pass
class DimensionalityError(Exception):
    """Raised in the event of dimensionality mismatch (e.g., three-dimensional mesh is requested for a two-dimensional geometry)."""
    pass
class DuplicateEntity(Exception):
    """Raised in the event that disallowed duplication is detected (e.g., blocks with the same name)."""
    pass

class Block:
    """An object that represents a contiguous structured grid.
    
    This object corresponds to the Plot3D notion of a block in that it defines either
    x,y or x,y,z coordinates in 2 or 3 dimensions in terms of 2 or 3 dimensional arrays.
    When initially created, it has no defined extents. Primitives are used to fill out
    the size of the arrays and which coordinates will be meshed, and grids are applied
    to determine the coordinates.

    Parameters
    ----------
    name:
        The name of the block.
    """
    def __init__(self, name:str):
        self.name = name
        self._primitives = []
        self.two_dimensional = False
        self._meshed = False
        self._indexed = False

    @classmethod
    def from_list(cls, name:str, array:list[Box], shape=None): # rework this with Numpy or something
        """Create a Block from a list of primitives.
        
        A class method that takes a list of primitives and creates a block from it. Currently only supports two dimensions.

        Parameters
        ----------
        name:
            The name of the block to create.
        array:
            The list of primitives to use.
        shape:
            The shape to use.  Passing None will cause the array to be treated as one-dimensional.

        Returns
        -------
        Block:
            The block created by arranging the primitives according to the shape.
        
        """
        primitives = []
        imax = len(array)
        jmax = 1
        # kmax = None
        if shape is not None:
            if len(shape) > 2:
                raise DimensionalityError("Three dimensional construction from array not allowed")
            imax = shape[0]
            jmax = shape[1]

        # Scan the rows to get the ni and nj values
        i_shift = [None] * imax
        j_shift = [None] * jmax
        index = 0
        for j in range(jmax):
            for i in range(imax):
                if array[index]:
                    if i_shift[i] is None:
                        i_shift[i] = array[index].ni
                    else:
                        if i_shift[i] != array[index].ni:
                            raise DimensionalityError('General dimensioning not allowed, all i sizes in array column must match')
                    if j_shift[j] is None:
                        j_shift[j] = array[index].nj
                    else:
                        if j_shift[j] != array[index].nj:
                            raise DimensionalityError('General dimensioning not implemented, all j sizes in array row must match')
                index += 1
        # Convert those into index shifts
        for i in range(1,imax-1):
            i_shift[i] += i_shift[i-1]
        i_shift.insert(0, 0)
        i_shift = i_shift[:-1]
        for j in range(1,jmax-1):
            j_shift[j] += j_shift[j-1]
        j_shift.insert(0, 0)
        j_shift = j_shift[:-1]
        # Now modify things so everything lines up
        index = 0
        for j in range(jmax):
            row = []
            for i in range(imax):
                row.append(array[index])
                #if array[index]:
                #    jmax = max(jmax, array[index].nj)
                index += 1
            for i,obj in enumerate(row):
                if not obj:
                    continue
                obj.i += i_shift[i]
                obj.j += j_shift[j]
                #if obj.nj != jmax:
                #    raise NotImplementedError
            #last_j += jmax
            primitives.extend([el for el in row if el is not None])
        object = cls(name)
        for primitive in primitives:
            object.add(primitive)

        return object

    @property
    def primitives(self):
        return self._primitives
    
    def add(self, primitive:Box):
        """Add a primitive to a Block.
        
        Add a primitive to a Block
        """
        if self._primitives:
            if primitive.two_dimensional == self.two_dimensional:
                self._primitives.append(primitive)
            else:
                return False
        else:
            self._primitives.append(primitive)
            self.two_dimensional = primitive.two_dimensional
        return True
    
    def mesh(self, xgrid, ygrid, zgrid=None): #, force=False):
        """Mesh a block with the given grids.
        
        Apply the x, y, and maybe z grid to a block, setting the coordinate locations.
        
        Parameters
        ----------
        xgrid:
            The x-direction grid.
        ygrid:
            The y-direction grid.
        zgrid: optional
            The z-direction grid, only required for three dimensional blocks.
            
        Raises
        ------
        AlreadyMeshed
            If the block is already meshed.
        """
        if self._meshed:
            raise AlreadyMeshed('Mesh "{self.name}" is already meshed')
        #if force:
        #    raise NotImplementedError
        
        x = np.zeros(self.ni+1)
        for i in range(self.ni+1):
            x[i] = xgrid.s(i)
        y = np.zeros(self.nj+1)
        for j in range(self.nj+1):
            y[j] = ygrid.s(j)

        if zgrid is not None:
            # Three dimensions
            if self.two_dimensional:
                raise DimensionalityError('z gridding not available for two-dimensional mesh')
            z = np.zeros(self.nk+1)
            for k in range(self.nk+1):
                z[k] = zgrid.s(k)
            
            self.x = np.zeros(self.node_count)
            self.y = np.zeros(self.node_count)
            self.z = np.zeros(self.node_count)

            k = 0
            index = 0
            for k in range(self.nk+1):
                for j in range(self.nj+1):
                    for i in range(self.ni+1):
                        if self.node_index[i, j, k] > 0:
                            self.x[index] = x[i]
                            self.y[index] = y[j]
                            self.z[index] = z[k]
                            index += 1
        else:
            # Two dimensions
            self.x = np.zeros(self.node_count)
            self.y = np.zeros(self.node_count)
            self.z = None

            k = 0
            index = 0
            for j in range(self.nj+1):
                for i in range(self.ni+1):
                    if self.node_index[i, j, k] > 0:
                        self.x[index] = x[i]
                        self.y[index] = y[j]
                        index += 1
        
        self._meshed = True

    def index(self): #, force=False):
        """Apply indices to the block based on primitives.
        
        Apply each of the primitives in mesh and determine the extents of the mesh.
        
        Raises
        ------
        AlreadyMeshed:
            If the block is already meshed.
        """
        if self._meshed:
            raise AlreadyMeshed('Mesh "{self.name}" is already meshed')
        #if force:
        #    raise NotImplementedError
        
        if not self._primitives:
            return
        
        # Get the blocks and sidesets we have
        self.blocks = []
        sidesets = set()
        for primitive in self._primitives:
            if primitive.element_set not in self.blocks:
                self.blocks.append(primitive.element_set)
            sidesets.update(primitive.sidesets())

        # Number the blocks from 1
        reverse_lookup = {}
        for i, block in enumerate(self.blocks):
            block.id = i+1
            reverse_lookup[block.id] = block
        
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

        self.cells = np.zeros((self.ni, self.nj, self.nk), dtype=np.uint8, order='F')
        self.cell_index = np.zeros((self.ni, self.nj, self.nk), dtype=np.uint64, order='F')
        node_nk = self.nk+1
        if self.two_dimensional:
            node_nk = 1
        self.node_index = np.zeros((self.ni+1, self.nj+1, node_nk), dtype=np.uint64, order='F')

        # Map it out
        for primitive in self.primitives:
            k0 = 0
            k1 = 1
            if not self.two_dimensional:
                k0 = primitive.k - self.k_offset
                k1 = primitive.k + primitive.nk - self.k_offset
            self.cells[primitive.i - self.i_offset:primitive.i + primitive.ni - self.i_offset,
                       primitive.j - self.j_offset:primitive.j + primitive.nj - self.j_offset,
                       k0:k1] = primitive.element_set.id * primitive.footprint()
        
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

        # Deal with sidesets
        self.sidesets = {}
        for set_name in sidesets:
            self.sidesets[set_name] = []
            for primitive in self._primitives:
                self.sidesets[set_name].extend(primitive.sideset(set_name))
                
        self._indexed = True
    
    def iblank(self) -> np.array:
        if not self._indexed:
            self.index()
        ib = np.where(self.node_index>0, 1, 0)
        return ib

    def write_exodusii(self, filename:str, info:list[str]|None=None):
        exo = exodusii.exodusii_file(filename, 'w')
        ndim = 3
        nnodes = 8
        type = 'HEX'
        if self.two_dimensional:
            ndim = 2
            nnodes = 4
            type = 'QUAD'
        exo.put_init(self.name, ndim, self.node_count, self.cell_count,
                     len(self.blocks), 0, len(self.sidesets))
        exo.put_coord(self.x, self.y, self.z)

        # Write out the blocks
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
            exo.put_element_conn(block.id, np.array(conn))

        # Write out the sidesets
        id = 1
        for name, subsets in self.sidesets.items():
            elements = []
            sides = []
            for sub in subsets:
                for k in range(sub.k - self.k_offset, sub.k - self.k_offset + sub.nk):
                    for j in range(sub.j - self.j_offset, sub.j - self.j_offset + sub.nj):
                        for i in range(sub.i - self.i_offset, sub.i - self.i_offset + sub.ni):
                            if self.cell_index[i,j,k] > 0:
                                elements.append(self.cell_index[i,j,k])
                                sides.append(sub.number)
            exo.put_side_set_param(id, len(elements))
            exo.put_side_set_name(id, name)
            exo.put_side_set_sides(id, elements, sides)
            id += 1

        if info is not None:
            # Assume that the lines are not too long
            exo.put_info(len(info), info)

        exo.close()
        #
        #exo = exodusii.exodusii_file(filename, 'r')
        #print(exo.get_info_records())
        #exo.close()
        return True
    
    def write_plot3d(self, filename:str, binary=True)->bool:
        # For now, only support one big Plot3D block
        if self.two_dimensional:
            the_shape = (self.ni+1, self.nj+1, 1)
            block = plot3d.Block(np.reshape(self.x, the_shape, order='F'),
                                 np.reshape(self.y, the_shape, order='F'),
                                 None)
            plot3d.write_plot2D(filename, [block], binary=binary)
        else:
            the_shape = (self.ni+1, self.nj+1, self.nk+1)
            block = plot3d.Block(np.reshape(self.x, the_shape, order='F'),
                                 np.reshape(self.y, the_shape, order='F'),
                                 np.reshape(self.z, the_shape, order='F'))
            plot3d.write_plot3D(filename, [block], binary=binary)
        return True
    
class Mesh:
    """An object that contains one or more grid blocks.
    
    This object contains the grid blocks and maintains the numbering and naming scheme.
    The other objects should be requested from this object, which enforces the necessary
    order of things and will minimize the potential issues with ids.

    Parameters
    ----------
    name:
        The name of the mesh.
    
    """
    def __init__(self, name:str):
        self.name = name
        self.blocks = []
        self.element_sets = {}

    def new_element_set(self, name:str)->ElementSet:
        """Create a new element set associated with this mesh.
        
        Parameters
        ----------
        name:
            Name of the new element set.

        Raises
        ------
        DupicateEntity:
            If the name is already present as an element set.

        """
        if name in self.element_sets:
            raise DuplicateEntity('An element set named "%s" is already part of mesh "%s", duplicates not allowed.' % (name, self.name))
        id = len(self.element_sets) + 1
        es = ElementSet(name, id)
        self.element_sets[name] = es
        return es
    
    def new_block(self, *primitives: Box):
        id = len(self.blocks) + 1
        block = Block('block_%d'%id)
        for primitive in primitives:
            if primitive.element_set is None:
                es_name = 'element_set_%d' % (len(self.element_sets) + 1)
                if primitive.name is not None:
                    es_name = 'element_set_' + primitive.name
                es = self.new_element_set(es_name)
                primitive.element_set = es
            elif isinstance(primitive.element_set, str):
                es = self.element_sets.get(primitive.element_set, self.new_element_set(primitive.element_set))
                primitive.element_set = es
            block.add(primitive)
        self.blocks.append(block)
        return block
    
    @property
    def node_count(self):
        count = 0
        for block in self.blocks:
            count += block.node_count
        return count
    
    @property
    def cell_count(self):
        count = 0
        for block in self.blocks:
            count += block.cell_count
        return count
    
    def block_from_list(self, name:str, array:list[Box], shape:tuple[int,int]|None=None)->Block:
        """Create a Block from a list of primitives.
        
        Take a list of primitives and create a block from it. Currently only supports two dimensions.

        Parameters
        ----------
        name:
            The name of the block to create.
        array:
            The list of primitives to use.
        shape:
            The shape to use.  Passing None will cause the array to be treated as one-dimensional.

        Returns
        -------
        Block:
            The block created by arranging the primitives according to the shape.
        
        """
        for primitive in array:
            if primitive is None:
                continue
            if primitive.element_set is None:
                es_name = 'element_set_%d' % (len(self.element_sets) + 1)
                if primitive.name is not None:
                    es_name = 'element_set_' + primitive.name
                primitive.element_set = self.new_element_set(es_name)
            elif isinstance(primitive.element_set, str):
                if primitive.element_set in self.element_sets:
                    primitive.element_set = self.element_sets[primitive.element_set]
                else:
                    primitive.element_set = self.new_element_set(primitive.element_set)
        block = Block.from_list(name, array, shape=shape)
        block.id = len(self.blocks) + 1
        self.blocks.append(block)
        return block
    
    def write_exodusii(self, filename:str, info:list[str]|None=None)->bool:
        return self.blocks[0].write_exodusii(filename, info=info)

    def write_plot3d(self, filename:str, binary=True)->bool:
        return self.blocks[0].write_plot3d(filename)
    
    def index(self):
        for block in self.blocks:
            # Handle the element sets
            for primitive in block.primitives:
                if primitive.element_set is None:
                    es_name = 'element_set_%d' % (len(self.element_sets) + 1)
                    if primitive.name is not None:
                        es_name = 'element_set_' + primitive.name
                    primitive.element_set = self.new_element_set(es_name)
                elif isinstance(primitive.element_set, str):
                    if primitive.element_set in self.element_sets:
                        primitive.element_set = self.element_sets[primitive.element_set]
                    else:
                        primitive.element_set = self.new_element_set(primitive.element_set)
            # Set up the block indices
            block.index()

    def mesh(self):
        for block in self.blocks:
            block.mesh()