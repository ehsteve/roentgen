X-ray emission from radionuclides
=================================

The purpose of this guide is to present an overview of the `roentgen.nuclides` module which provides access to tables of emission lines from radioactive sources (e.g. Am-241, Ba-133).
radionuclides are unstable atoms that spontaneously and randomly decay into a different nuclide over time, emitting ionizing radiation in the process.
Those new nuclide may undergo a similar process leading to a chain of decay of different elements.
The radiation emitted may taken the form of x-rays or gamma-rays with characteristic energies.
Because these characteristic energies are well-known, these sources are frequently used in laboratories to calibrate x-ray and gamma-ray detectors.

This module provides access to data from the Library for gamma and alpha emissions (`Lara <http://www.lnhb.fr/accueil/donnees-nucleaires/module-lara/>`_) provided by the `Laboratoire National Henri Becquerel <http://www.lnhb.fr>`_.

The list of available nuclides is available through::

    >>> from roentgen.nuclides import nuclides_list
    >>> nuclides_list
    <QTable length=351>
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
    >>> print(fe55.lines)
    <QTable length=5>
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
    >>> print(strong_lines)
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

The Nuclide class also provides access to some helpful information such as the half life

    >>> print(am241.half_life)
    <Quantity 432.60577484 yr>

and the url for the data sheet

    >>> print(am241.data_sheet)
    'http://www.lnhb.fr/nuclides/Am-241_tables.pdf'
