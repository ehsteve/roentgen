"""
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show app.py
in your browser.
"""
import numpy as np
from roentgen.absorption import Material, get_density
import astropy.units as u
import roentgen

from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models.widgets import (
    Slider,
    TextInput,
    Select,
    Paragraph,
    AutocompleteInput,
    CheckboxGroup,
    TableColumn,
    DataTable,
)
from bokeh.plotting import figure

R = 287.058 * u.J / u.kg / u.Kelvin
from astropy import constants as const


def get_air_density(pressure, temperature):
    return pressure / (R * temperature.to('K', equivalencies=u.temperature()))

DEFAULT_MATERIAL = ["Si", "Si", "Si", "Si", "Si"]
DEFAULT_THICKNESS = [100.0, 0, 0, 0, 0] * u.micron
DEFAULT_ENERGY_LOW = 1.0
DEFAULT_ENERGY_HIGH = 50.0
DEFAULT_DENSITY = [2.33, 2.33, 2.33, 2.33, 2.33]
DEFAULT_TYPE = "Transmission"
NUMBER_OF_MATERIALS = len(DEFAULT_MATERIAL)

DEFAULT_DETECTOR_MATERIAL = 'cdte'
DEFAULT_DETECTOR_THICKNESS = 1 * u.mm
DEFAULT_DETECTOR_DENSITY = get_density(DEFAULT_DETECTOR_MATERIAL)

DEFAULT_AIR_THICKNESS = 1 * u.m
DEFAULT_AIR_PRESSURE = 1 * const.atm
DEFAULT_AIR_TEMPERATURE = 25 * u.Celsius

PLOT_HEIGHT = 300
PLOT_WIDTH = 900
TOOLS = "pan,box_zoom,box_select,crosshair,undo,redo,save,reset,hover"

# defaults
material_list = []
energy = u.Quantity(np.arange(DEFAULT_ENERGY_LOW, DEFAULT_ENERGY_HIGH, 1), "keV")

this_material = Material(DEFAULT_MATERIAL[0], DEFAULT_THICKNESS[0])
y = this_material.transmission(energy).value / 100.0
air_density = get_air_density(DEFAULT_AIR_PRESSURE, DEFAULT_AIR_TEMPERATURE)
air = Material('air', DEFAULT_AIR_THICKNESS, density=air_density)
y *= air.transmission(energy).value / 100.0
this_detector = Material(DEFAULT_DETECTOR_MATERIAL, DEFAULT_DETECTOR_THICKNESS)
y *= this_detector.absorption(energy).value / 100.0

x = energy.value
source = ColumnDataSource(data={"x": x, "y": y * 100})

all_materials = list(roentgen.elements["symbol"]) + \
                list(roentgen.compounds["symbol"])
all_materials.sort()
all_materials = [this_material.lower() for this_material in all_materials]

# Set up the plot
plot = figure(
    plot_height=PLOT_HEIGHT,
    plot_width=PLOT_WIDTH,
    title=DEFAULT_TYPE,
    tools=TOOLS,
    x_range=[DEFAULT_ENERGY_LOW, DEFAULT_ENERGY_HIGH],
    y_range=[0, 100],
)
plot.yaxis.axis_label = "%"
plot.xaxis.axis_label = "Energy [keV]"
plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)

plot.title.text = "boo"

# Set up the inputs
ylog_checkbox = CheckboxGroup(labels=["y-log"], active=[0])

# Set plot widgets
energy_low = Slider(title="energy (low)", value=3, start=1, end=5.0, step=0.1)
energy_high = Slider(title="energy (high)", value=50, start=10, end=100, step=1)
energy_step = Slider(title="energy steps", value=1, start=0.1, end=2, step=0.1)

# materials in the path
material_input = TextInput(title="Material", value=DEFAULT_MATERIAL[0])
thickness_input = TextInput(title="thickness", value=f"{DEFAULT_THICKNESS[0]}")
density_input = TextInput(title="density", value=f"{this_material.density}")

air_thickness_input = TextInput(title="air path length", value=f"{DEFAULT_AIR_THICKNESS}")
air_pressure_input = TextInput(title='air pressure', value=f"{DEFAULT_AIR_PRESSURE}")
air_temperature_input = TextInput(title='air temperature', value=f"{DEFAULT_AIR_TEMPERATURE}")

detector_material_input = TextInput(title='Detector', value=DEFAULT_DETECTOR_MATERIAL)
detector_thickness_input = TextInput(title='thickness', value=str(DEFAULT_DETECTOR_THICKNESS))
detector_density_input = TextInput(title='density', value=str(DEFAULT_DETECTOR_DENSITY))

p = Paragraph(text="Nothing to see here", width=200)
p.text = "All available materials: " + ", ".join(all_materials)

columns = [
    TableColumn(field="x", title="energy [keV]"),
    TableColumn(field="y", title="Percent"),
]
data_table = DataTable(source=source, columns=columns, width=400, height=700)


def update_data(attrname, old, new):

    # if ylog_checkbox.active:
    #    plot.y_axis_type = "log"
    # else:
    #    plot.y_axis_type = "auto"

    energy = u.Quantity(
        np.arange(float(energy_low.value), float(energy_high.value),
                  float(energy_step.value)), "keV")
    x = energy.value
    i = 0

    material_list = material_input.value.split(',')
    print(material_list)
    plot_title = ""

    for this_material_str in material_list:
        this_material = Material(str(this_material_str),
                                 u.Quantity(thickness_input.value).to("micron"),
                                 density=u.Quantity(density_input.value).to("g / cm ** 3"))
        if i == 0:
            y = this_material.transmission(energy).value / 100.0
        else:
            y *= this_material.transmission(energy).value / 100.0
        i += 0
        plot_title += f"{this_material.name} {this_material.thickness}"

    # parse atm input
    if air_pressure_input.value.count('atm') > 0:
        air_pressure = air_pressure_input.value.split('atm')[0] * const.atm
    else:
        air_pressure = u.Quantity(air_pressure_input.value)
    air_density = get_air_density(air_pressure,
                                  u.Quantity(air_temperature_input.value))
    air = Material('air', u.Quantity(air_thickness_input.value).to("meter"),
                   density=air_density)
    y *= air.transmission(energy).value / 100.0
    plot_title += f" {air.name} {air.density.to('g / cm**3'):.2E} {air.thickness:.2f}"

    this_detector = Material(detector_material_input.value, u.Quantity(detector_thickness_input.value),
                             density=u.Quantity(detector_density_input.value))
    y *= this_detector.absorption(energy).value / 100.0
    plot_title += f" {this_detector.name} {this_detector.thickness:.2f}"
    plot.title.text = plot_title
    source.data = dict(x=x, y=y * 100)


def update_detector_density_and_data(attrname, old, new):
    detector_density_input.value = str(get_density(detector_material_input.value))
    update_data(attrname, old, new)
    return attrname, old, new

# plot.x_range.on_change('start', update_data)
# plot.x_range.on_change('end', update_data)

update_input_list = [material_input, thickness_input, density_input, air_thickness_input,
                     air_pressure_input, air_temperature_input, detector_density_input,
                     detector_thickness_input, detector_material_input, energy_low,
                     energy_step, energy_high]

for w in update_input_list:
    w.on_change("value", update_data)

detector_material_input.on_change("value", update_detector_density_and_data)
#material_input.on_change("value", update_material_fields)

density_input.disabled = True
detector_density_input.disabled = True

curdoc().add_root(
    layout(
        [
            [material_input, thickness_input, density_input],
            [air_pressure_input, air_thickness_input, air_temperature_input],
            [detector_material_input, detector_thickness_input, detector_density_input],
            # [material_input_list[3], thickness_input_list[3], density_input_list[3]],
            # [material_input_list[4], thickness_input_list[4], density_input_list[4]],
            [energy_low, energy_high, energy_step],
            [plot, data_table, p],
        ],
        sizing_mode="fixed",
    )
)
curdoc().title = "Roentgen"
