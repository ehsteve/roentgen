import astropy.units as u

import roentgen.util as util
from roentgen.absorption import Material

m = Material({"Si": 1, "He": 1}, 5 * u.m)

stack1 = Material("Ge", 500 * u.micron) + Material("cdte", 100 * u.micron)
stack2 = Material("Si", 500 * u.micron) + Material("Al", 500 * u.micron)

total_stack = stack1 + stack2

print(m.density)
d = util.get_material_density("Si") * 0.5 + util.get_material_density("He") * 0.5
print(d == m.density)
