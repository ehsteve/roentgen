# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.
import os
import pytest
import astropy.units as u

packagename = os.path.basename(os.path.dirname(__file__))
