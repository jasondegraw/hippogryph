# SPDX-FileCopyrightText: 2014 Jason W. DeGraw <jason.degraw@gmail.com>
# SPDX-FileCopyrightText: 2023-present Oak Ridge National Laboratory, managed by UT-Battelle
#
# SPDX-License-Identifier: BSD-3-Clause
import hippogryph

def dev_null(mesg):
    pass

def test_one_sided_vinokur_functions():
    ds = hippogryph.single_sided_vinokur(1.0e-6, 1, 16, output=dev_null)
    assert abs(ds - 7.4717023123662765) < 1.0e-15
    assert abs(hippogryph.vkruh(ds, 1.0, 1, 16) - 1.0e-6) < 1.0e-15

    ds = hippogryph.single_sided_vinokur(-1.0e-6, 1, 16, output=dev_null)
    assert ds is None

    ds = hippogryph.single_sided_vinokur(2.0, 1, 16, output=dev_null)
    assert ds is None

    ds = hippogryph.single_sided_vinokur(1.0e-6, 1, 16, output=dev_null, max_iterations=2)
    assert ds is None

def test_one_sided_vinokur_object():
    vkr = hippogryph.VinokurSingleSided.from_delta(1.0e-6, 1.0, 1, 16, output=dev_null)
    assert abs(vkr.s(1) - 1.0e-6) < 1.0e-15
    assert abs(vkr.s(16) - 1.0) < 1.0e-15

    vkr = hippogryph.VinokurSingleSided.from_delta(2.0, 1.0, 1, 16, output=dev_null)
    assert vkr is None

def test_uniform_object():
    obj = hippogryph.Uniform.from_delta(0.1, 10)
    assert obj.L == 1.0
    assert obj.s(1) == obj.delta
    obj = hippogryph.Uniform.from_intervals(1.0, 10)
    assert obj.delta == 0.1
    assert obj.s(1) == obj.delta

def test_geometric_function():
    factor = hippogryph.single_sided_geometric(1.0e-3, 32, output=dev_null)
    assert abs(hippogryph.geometric(factor, 1.0e-3, 32) - 1.0) < 1.0e-14
    assert abs(hippogryph.geometric_sum(factor, 1.0e-3, 32) - 1.0) < 1.0e-14
    assert hippogryph.geometric(factor, 1.0e-3, 0) == 0.0
    assert hippogryph.geometric(factor, 1.0e-3, 1) == 1.0e-3
    assert hippogryph.geometric(1.0, 0.1, 2) == 0.2
    assert hippogryph.geometric_sum(factor, 1.0e-3, 0) == 0.0
    assert hippogryph.geometric_sum(factor, 1.0e-3, 1) == 1.0e-3

def test_geometric_object():
    obj = hippogryph.Geometric.from_delta(0.025, 10.0, 32, output=dev_null)
    #assert abs(obj.factor - 11.328501284496937) < 1.0e-14
    assert obj.s(1) == 0.025
    assert obj.s(2) == 0.025 * (1.0 + obj.factor)
    assert abs(obj.s(32) - 10.0) < 1.0e-14

def test_composites():
    # Two uniform grids
    grids = [hippogryph.Uniform.from_delta(0.1, 10),
             hippogryph.Uniform.from_delta(0.2, 5)]
    obj = hippogryph.Composite(grids)
    assert obj.L == 2.0
    assert obj.N == 15
    assert obj.s(10) == 1.0
    assert obj.s(11) == 1.2

    grids[1] = hippogryph.Geometric.from_delta(0.025, 10.0, 32, output=dev_null)
    obj = hippogryph.Composite(grids)
    assert obj.L == 11.0
    assert obj.N == 42
    assert obj.s(42) == 11.0
    assert obj.intervals == [10, 42]

    grids = []
    obj = hippogryph.Composite(grids)
    assert obj.grids == []

def test_fails():
    factor = hippogryph.single_sided_geometric(-0.002, 32, output=dev_null)
    assert factor is None
    factor = hippogryph.single_sided_geometric(2.0, 32, output=dev_null)
    assert factor is None
    grid = hippogryph.Geometric.from_delta(25.0, 10.0, 0, 32, output=dev_null)
    assert grid is None
    factor = hippogryph.single_sided_geometric(1.0e-3, 32, output=dev_null, max_iterations=2)
    assert factor is None