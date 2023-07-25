# SPDX-FileCopyrightText: 2023-present Jason W. DeGraw <jason.degraw@gmail.com>
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