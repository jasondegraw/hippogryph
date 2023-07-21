# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause
import argparse
import hippogryph

def main(inputs = None):
    if inputs is None:
        inputs = {}
    factor = hippogryph.single_sided_factor(1.0e-6, 1, 16)
    print(factor)
    for i in range(1, 17):
        print(i, hippogryph.vkruh(factor, 1.0, i, 16))