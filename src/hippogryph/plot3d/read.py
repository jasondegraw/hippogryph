# SPDX-License-Identifier: 	NASA-1.3

import numpy as np 
import os.path as osp
import struct
from typing import List
from .block import Block

def __read_plot3D_chunk_binary(f,IMAX:int,JMAX:int,KMAX:int, big_endian:bool=False,read_double:bool=True):
    """Reads and formats a binary chunk of data into a plot3D block

    Args:
        f (io): file handle
        IMAX (int): maximum I index
        JMAX (int): maximum J index
        KMAX (int): maximum K index
        big_endian (bool, Optional): Use big endian format for reading binary files. Defaults False.

    Returns:
        numpy.ndarray: Plot3D variable either X,Y, or Z 
    """
    A = np.empty(shape=(IMAX, JMAX, KMAX))
    for k in range(KMAX):
        for j in range(JMAX):
            for i in range(IMAX):
                if read_double:
                    A[i,j,k] = struct.unpack(">d",f.read(8))[0] if big_endian else struct.unpack("<d",f.read(8))[0]
                else:
                    A[i,j,k] = struct.unpack(">f",f.read(4))[0] if big_endian else struct.unpack("<f",f.read(4))[0]
    return A

def read_word(f):
    """Continously read a word from an ascii file

    Args:
        f (io): file handle

    Yields:
        float: value from ascii file
    """
    for line in f:
        line = line.strip().replace('\n','').split(' ')
        tokenArray = [float(entry) for entry in line if entry]
        for token in tokenArray:
            yield token

def __read_plot3D_chunk_ASCII(f,IMAX:int,JMAX:int,KMAX:int):
    """Reads and formats a binary chunk of data into a plot3D block

    Args:
        f (io): file handle
        IMAX (int): maximum I index
        JMAX (int): maximum J index
        KMAX (int): maximum K index
        big_endian (bool, Optional): Use big endian format for reading binary files. Defaults False.

    Returns:
        numpy.ndarray: Plot3D variable either X,Y, or Z 
    """
    tokenArray = np.zeros(shape=(IMAX*JMAX*KMAX))
    i = 0
    for w in read_word(f):
        tokenArray[i] = w
        i+=1
        if i>len(tokenArray)-1:
            break

    A = np.reshape(tokenArray,newshape=(KMAX,JMAX,IMAX))
    A = np.transpose(A,[2,1,0])    
    return A

def read_plot3D(filename:str, binary:bool=True,big_endian:bool=False,read_double:bool=True):
    """Reads a plot3d file and returns Blocks

    Args:
        filename (str): name of the file to read, .p3d, .xyz, .pdc, .plot3d? 
        binary (bool, optional): indicates if the file is binary. Defaults to True.
        big_endian (bool, optional): use big endian format for reading binary files
        read_float (bool, optional): read floating point. Only affects binary files

    Returns:
        List[Block]: List of blocks insdie the plot3d file
    """
    
    blocks = list()
    if osp.isfile(filename):
        if binary:
            with open(filename,'rb') as f:
                nblocks = struct.unpack(">I",f.read(4))[0] if big_endian else struct.unpack("I",f.read(4))[0] # Read bytes            
                IMAX = list(); JMAX = list(); KMAX = list()
                for b in range(nblocks):
                    if big_endian:
                        IMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        KMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                    else:
                        IMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                        KMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes

                for b in range(nblocks):
                    X = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Y = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Z = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    b_temp = Block(X,Y,Z)                    
                    blocks.append(b_temp)
        else:
            with open(filename,'r') as f: 
                nblocks = int(f.readline())
                IMAX = list(); JMAX = list(); KMAX = list()
                
                for b in range(nblocks):
                    IJK = f.readline().replace('\n','').split(' ')
                    tokens = [int(w) for w in IJK if w]
                    IMAX.append(tokens[0])
                    JMAX.append(tokens[1])
                    KMAX.append(tokens[2])            

                for b in range(nblocks):
                    X = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Y = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Z = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    b_temp = Block(X,Y,Z)                    
                    blocks.append(b_temp)
    return blocks

def read_plot2D(filename:str, binary:bool=True,big_endian:bool=False,read_double:bool=True):
    """Reads a plot2d file and returns Blocks

    Args:
        filename (str): name of the file to read, .p3d, .xyz, .pdc, .plot3d? 
        binary (bool, optional): indicates if the file is binary. Defaults to True.
        big_endian (bool, optional): use big endian format for reading binary files
        read_float (bool, optional): read floating point. Only affects binary files

    Returns:
        List[Block]: List of blocks insdie the plot3d file
    """
    
    blocks = list()
    if osp.isfile(filename):
        if binary:
            with open(filename,'rb') as f:
                nblocks = struct.unpack(">I",f.read(4))[0] if big_endian else struct.unpack("I",f.read(4))[0] # Read bytes            
                IMAX = list(); JMAX = list(); KMAX = list()
                for b in range(nblocks):
                    if big_endian:
                        IMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                    else:
                        IMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                    KMAX.append(1)

                for b in range(nblocks):
                    X = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Y = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Z = None
                    b_temp = Block(X,Y,Z)                    
                    blocks.append(b_temp)
        else:
            with open(filename,'r') as f: 
                nblocks = int(f.readline())
                IMAX = list(); JMAX = list(); KMAX = list()
                
                for b in range(nblocks):
                    IJK = f.readline().replace('\n','').split(' ')
                    tokens = [int(w) for w in IJK if w]
                    IMAX.append(tokens[0])
                    JMAX.append(tokens[1])
                    KMAX.append(1)

                for b in range(nblocks):
                    X = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Y = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Z = None
                    b_temp = Block(X,Y,Z)                    
                    blocks.append(b_temp)
    return blocks


