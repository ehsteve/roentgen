import pytest

from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.material import MassAttenuationCoefficient


def test_interpolate_matches_at_data():
    for this_element in ['Te', 'Si', 'Ge', 'cdte']:
        te = MassAttenuationCoefficient(this_element)
        assert_allclose(te.data, te.func(te.energy))

        # now change one point and make sure it no longer works
        te.energy[0] = 10 * u.keV
        with pytest.raises(AssertionError):
            assert_allclose(te.data, te.func(te.energy))
