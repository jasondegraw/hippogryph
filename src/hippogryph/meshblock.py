# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

class Block:
    def __init__(self, name, i, j, k=None):
        self.name = name
        self.i = i
        self.j = j
        self.k = k
        self.two_dimensional = k is None

class Mesh:
    def __init__(self):
        self._blocks = []
        self.two_dimensional = False

    @property
    def blocks(self):
        return self._blocks
    
    def add(self, block):
        if self._blocks:
            if block.two_dimensional == self.two_dimensional:
                self._blocks.append(block)
            else:
                return False
        else:
            self._blocks.append(block)
            self.two_dimensional = block.two_dimensional
