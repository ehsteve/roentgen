X-ray emission from radionuclides
=================================

.. warning::

    This module is in beta and its API may change in future releases.
    It does not yet calculate the full decay chain or account for possible emissions from decay products.
    The emission from radionuclides which decay straight into stable nuclides should be accurate.

The purpose of this guide is to present an overview of the `roentgen.nuclides` module which provides access to tables of emission lines from radioactive sources (e.g. Am-241, Ba-133).
Radionuclides are unstable atoms that spontaneously and randomly decay into a different nuclide over time, emitting ionizing radiation in the process.
Those new nuclide may undergo a similar process leading to a chain of decay of different elements.
The radiation emitted may take the form of x-rays or gamma-rays with characteristic energies.
Because these characteristic energies are well-known and sealed sources of radiation are readily available, these sources are frequently used in laboratories to calibrate x-ray and gamma-ray detectors.

This module provides access to data from the Library for gamma and alpha emissions (`Lara <http://www.lnhb.fr/accueil/donnees-nucleaires/module-lara/>`_) provided by the `Laboratoire National Henri Becquerel <http://www.lnhb.fr>`_.

The list of available nuclides is available through::

    >>> from roentgen.nuclides import nuclides_list
    >>> print(nuclides_list) # doctest: +SKIP
    symbol     name    mass_number       filename            half_life
                                                                yr
    str2     str12       int64           str19               float64
    ------ ----------- ----------- ------------------- ----------------------
        H    Hydrogen           3    H-3_@00.lara.txt     12.310822115750247
        He      Helium           6   He-6_@00.lara.txt 2.5575772555580905e-08
        Be   Beryllium           7   Be-7_@00.lara.txt     0.1457018277689051
        Be   Beryllium          10  Be-10_@00.lara.txt     1393008.3403047125
        C      Carbon          11   C-11_@00.lara.txt  3.871333688239917e-05
        N    Nitrogen          13   N-13_@00.lara.txt  1.895011027454559e-05
    ...         ...         ...                 ...                    ...
        Am   Americium         244 Am-244_@00.lara.txt  0.0011521788729180924
        Pu   Plutonium         244 Pu-244_@00.lara.txt      81089816.71610008
        Cm      Curium         245 Cm-245_@00.lara.txt      8248.409257991734
        Am   Americium         245 Am-245_@00.lara.txt 0.00023385808806753364
        Cm      Curium         246 Cm-246_@00.lara.txt      4721.525084290313
        Cm      Curium         248 Cm-248_@00.lara.txt     347935.20419803786
        Cf Californium         252 Cf-252_@00.lara.txt      2.646905975105838

To download the entire list :download:`nuclides_list.csv <../../roentgen/data/nuclides_list.csv>` or view it online

To access a particular nuclide first initialize it::

    >>> from roentgen.nuclides import Nuclide
    >>> fe55 = Nuclide('Fe', 55)
    >>> print(fe55.lines) # doctest: +SKIP
    energy intensity  origin parent
    keV
    float64  float64    str7   str7
    ------- --------- ------- ------
    0.6635      0.523     Mn   Fe-55
    5.88772       8.4     Mn   Fe-55
    5.89881     16.48     Mn   Fe-55
    6.513        3.38     Mn   Fe-55
    125.949   1.3e-07  Mn-55   Fe-55

For nuclides with complicated decay chains this may be too much information.
This module therefore provides a convenience function to filter the list.
For example, Am-241 has a long decay chain generating 873 different lines some of them very weak.
We can there filter the list to only the strongest lines within a particular energy range:::

    >>> import astropy.units as u
    >>> am241 = Nuclide('Am', 241)
    >>> strong_lines = am241.get_lines(5 * u.keV, 100 * u.keV, min_intensity=5)
    >>> print(strong_lines) # doctest: +SKIP
    energy  intensity  origin parent
    keV
    -------- --------- ------- ------
    12.20145      8.04     Pb  Tl-209
    13.852       13.02     Np  Am-241
    14.08955      19.0     Fr  Ac-225
    14.89645      13.6     Ac  Ra-225
    15.7405       59.7     Pa  Np-237
    16.1665       40.6      U  Pa-233
    16.96        18.58     Np  Am-241
    29.374        14.3  Pa-233 Np-237
    40.09         30.0  Ac-225 Ra-225
    59.5409      35.92  Np-237 Am-241
    72.8049       5.85     Pb  Tl-209
    74.97         9.84     Pb  Tl-209
    86.477       12.26  Pa-233 Np-237
    94.666         9.1      U  Pa-233
    98.44         14.6      U  Pa-233

The Nuclide class also provides access to some additional information such as the half life

    >>> print(am241.half_life)
    432.6057748371232 yr

and the url for the data sheet

    >>> print(am241.data_sheet)
    http://www.lnhb.fr/nuclides/Am-241_tables.pdf

All other metadata is available in meta attribute

    >>> print(am241.meta)
    {'Nuclide': 'Am-241', 'Element': 'Americium', 'Z': 95.0, 'Daughter(s)': [{'decay_mode': 'alpha', 'element': 'Np', 'mass_number': 237, 'branching_ratio': 100.0}], 'Qalpha': 5637.82, 'Possible parent(s)': '(B-)', 'Half-life (a)': <Quantity 432.6 yr>, 'Half-life (s)': <Quantity 1.3652e+10 s>, 'Decay constant (1/s)': <Quantity 5.077e-11 1 / s>, 'Specific activity (Bq/g)': <Quantity 1.2688e+11 Bq / g>, 'Reference': 'KRI - 2009'}

Ba-133 Example Spectrum
------------------------

As an example, let's consider the nuclide Ba-133. This nuclide decays through electron capture into Cs-133, emitting a number of gamma-rays and x-rays in the process.
The strongest emissions are at 31 keV, 81 keV, and 356 keV.

    >>> from roentgen.nuclides import Nuclide
    >>> ba133 = Nuclide('Ba', 133)
    >>> print(ba133)
    Nuclide: Ba-133, (Barium) half life=10.538824245189748 yr - (14 lines)
    Daughters: [{'decay_mode': 'EC', 'element': 'Cs', 'mass_number': 133, 'branching_ratio': 100.0}]
    <QTable length=14>
    energy  intensity  origin parent
    keV
    float64   float64    str7   str7
    -------- --------- ------- ------
    4.67355     15.87  Cs-133 Ba-133
    30.6254      33.8  Cs-133 Ba-133
    30.9731      62.4  Cs-133 Ba-133
    35.053     18.24  Cs-133 Ba-133
    35.9003      4.45  Cs-133 Ba-133
    53.1622      2.14  Cs-133 Ba-133
    79.6142      2.63  Cs-133 Ba-133
    80.9979     33.31  Cs-133 Ba-133
    160.6121     0.638  Cs-133 Ba-133
    223.2368      0.45  Cs-133 Ba-133
    276.3989      7.13  Cs-133 Ba-133
    302.8508     18.31  Cs-133 Ba-133
    356.0129     62.05  Cs-133 Ba-133
    383.8485      8.94  Cs-133 Ba-133

With this information, it is possible to simulate the emission from a radionuclide as it would be seen by a detector with a finite energy resolution.
In the following example, we will simulate the emission from a Ba-133 source as it would be seen by a detector with a 0.5 keV energy resolution up to 80 keV.
Given this energy resolution, the individual lines will be broadened and will overlap.

.. plot::
    :include-source:

    import numpy as np

    import astropy.units as u
    from astropy.modeling.models import Gaussian1D
    from astropy import visualization

    from matplotlib import pyplot as plt
    from roentgen.nuclides import Nuclide

    ba133 = Nuclide('Ba', 133)
    lines = ba133.get_lines(1 * u.keV, 100 * u.keV, min_intensity=1)
    energy_resolution = 0.5 * u.keV

    g0 = Gaussian1D(amplitude=lines[0]['intensity'], mean=lines[0]['energy'], stddev=energy_resolution)
    for this_line in lines[1:]:
        g0 += Gaussian1D(amplitude=this_line['intensity'], mean=this_line['energy'], stddev=energy_resolution)

    energy_ax = np.arange(1, 85, 0.5) * u.keV

    with visualization.quantity_support():
        fig, ax = plt.subplots()
        ax.plot(energy_ax, g0(energy_ax))
        ax.set_xlabel(f'Energy')
        ax.vlines(lines['energy'], ymin=[0]*len(lines), ymax=lines['intensity'])
        plt.show()