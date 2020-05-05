import astropy.units as u

import roentgen

__all__ = ["is_an_element",
           "get_atomic_number",
           "is_in_known_compounds",
           "get_compound_index",
           "get_density",
           "get_element_symbol",
           "density_ideal_gas"]


def is_an_element(element_str):
    """Returns True if the string represents an element"""
    result = False
    lower_case_list = list([s.lower() for s in roentgen.elements["symbol"]])
    if (len(element_str) <= 2) and (element_str.lower() in lower_case_list):
        result = True
    else:
        lower_case_list = list([s.lower() for s in roentgen.elements["name"]])
        if element_str.lower() in lower_case_list:
            result = True
    return result


def get_element_symbol(element_str):
    """Return the element abbreviation"""
    lower_case_symbol_list = list([s.lower() for s in roentgen.elements["symbol"]])
    if element_str.lower() in lower_case_symbol_list:  # already a symbol
        return element_str.capitalize()
    elif is_an_element(element_str):
        lower_case_list = list([s.lower() for s in roentgen.elements["name"]])
        return roentgen.elements[lower_case_list.index(element_str.lower())]["symbol"].capitalize()
    else:
        return None


def get_atomic_number(element_str):
    """Return the atomic number of the element"""
    # check to see if element_str is symbol
    if is_an_element(element_str):
        if len(element_str) <= 2:
            lower_case_list = list([s.lower() for s in roentgen.elements["symbol"]])
            atomic_number = roentgen.elements[
                lower_case_list.index(element_str.lower())
            ]["z"]
        else:
            lower_case_list = list([s.lower() for s in roentgen.elements["name"]])
            atomic_number = roentgen.elements[
                lower_case_list.index(element_str.lower())
            ]["z"]
    else:
        atomic_number = None
    return atomic_number


def is_in_known_compounds(compound_str):
    """Returns True is the compound is in the list of known compounds"""
    lcase_symbols_list = list([s.lower() for s in roentgen.compounds["symbol"]])
    lcase_name_list = list([s.lower() for s in roentgen.compounds["name"]])
    case1 = compound_str.lower() in lcase_symbols_list
    case2 = compound_str.lower() in lcase_name_list
    return case1 or case2


def get_compound_index(compound_str):
    """Return the index of the compound in the compound table"""
    if is_in_known_compounds(compound_str):
        lower_case_symbols_list = list(
            [s.lower() for s in roentgen.compounds["symbol"]]
        )
        if compound_str.lower() in lower_case_symbols_list:
            return lower_case_symbols_list.index(compound_str.lower())
        lcase_name_list = list([s.lower() for s in roentgen.compounds["name"]])
        if compound_str.lower() in lcase_name_list:
            return lcase_name_list.index(compound_str.lower())
    else:
        return None


def get_density(material_str):
    """Given a material name return the default density"""
    if is_an_element(material_str):
        ind = get_atomic_number(material_str) - 1
        density = roentgen.elements[ind]["density"]
    else:
        # not using loc because table indexing is not yet stable
        # self.density = roentgen.compounds.loc[material_str]['density']
        index = list(roentgen.compounds["symbol"]).index(material_str)
        density = roentgen.compounds[index]["density"]
    return density


u.quantity_input(pressure=u.pascal)
def density_ideal_gas(pressure, temperature):  # noqa
    """Given pressure and temperature of a dry gas, return the density using
    the ideal gas law"""
    R = 287.058 * u.J / u.kg / u.Kelvin
    return pressure / (R * temperature.to('K', equivalencies=u.temperature()))
