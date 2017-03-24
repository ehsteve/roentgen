from roentgen.absorption import Material
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

thickness = 500 * u.micron

material_list = ['cdte', 'Si', 'Ge']
energy = u.Quantity(np.arange(1, 100, 0.1), 'keV')

for material in material_list:
    mat = Material(material, thickness)
    plt.plot(energy, mat.absorption(energy), label=mat.name)

plt.xlabel('Energy [' + str(energy.unit) + ']')
plt.ylabel('Efficiency')
plt.legend(loc='lower left')
plt.show()
