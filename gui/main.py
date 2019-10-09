"""
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show gui
in your browser.
"""
import numpy as np
from os.path import dirname, join

from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource, CustomJS
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
    Button,
    CheckboxGroup,
    Toggle,
    Div
)
from bokeh.plotting import figure

import astropy.units as u
from astropy import constants as const

from roentgen.absorption import Material, get_density
import roentgen

R = 287.058 * u.J / u.kg / u.Kelvin


def get_air_density(pressure, temperature):
    return pressure / (R * temperature.to('K', equivalencies=u.temperature()))

DEFAULT_MATERIAL = ["Si", "Si", "Si", "Si", "Si"]
DEFAULT_THICKNESS = [100.0, 0, 0, 0, 0] * u.micron
DEFAULT_ENERGY_LOW = 1.0
DEFAULT_ENERGY_HIGH = 200.0
DEFAULT_DENSITY = [2.33, 2.33, 2.33, 2.33, 2.33]
DEFAULT_TYPE = "Transmission"
NUMBER_OF_MATERIALS = len(DEFAULT_MATERIAL)

DEFAULT_DETECTOR_MATERIAL = 'cdte'
DEFAULT_DETECTOR_THICKNESS = 1 * u.mm
DEFAULT_DETECTOR_DENSITY = get_density(DEFAULT_DETECTOR_MATERIAL)

DEFAULT_ENERGY_RESOLUTION = 0.25

DEFAULT_AIR_THICKNESS = 1 * u.m
DEFAULT_AIR_PRESSURE = 1 * const.atm
DEFAULT_AIR_TEMPERATURE = 25 * u.Celsius

PLOT_HEIGHT = 300
PLOT_WIDTH = 900
TOOLS = "pan,box_zoom,box_select,crosshair,undo,redo,save,reset,hover"

# defaults
material_list = []
energy = u.Quantity(np.arange(DEFAULT_ENERGY_LOW, DEFAULT_ENERGY_HIGH,
                              DEFAULT_ENERGY_RESOLUTION), "keV")

this_material = Material(DEFAULT_MATERIAL[0], DEFAULT_THICKNESS[0])
y = this_material.transmission(energy)
air_density = get_air_density(DEFAULT_AIR_PRESSURE, DEFAULT_AIR_TEMPERATURE)
air = Material('air', DEFAULT_AIR_THICKNESS, density=air_density)
y *= air.transmission(energy)
this_detector = Material(DEFAULT_DETECTOR_MATERIAL, DEFAULT_DETECTOR_THICKNESS)
y *= this_detector.absorption(energy)

x = energy.value
source = ColumnDataSource(data={"x": x, "y": y})

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
    x_range=[1, 50],
    y_range=[0, 1],
)
plot.yaxis.axis_label = "fraction"
plot.xaxis.axis_label = "Energy [keV]"
plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)
plot.title.text = ""

# Set up the inputs
ylog_checkbox = CheckboxGroup(labels=["y-log"], active=[0])

# Set plot widgets
#energy_low = Slider(title="energy (low) [keV]", value=DEFAULT_ENERGY_LOW,
#                    start=1, end=5.0, step=0.1)
#energy_high = Slider(title="energy (high) [keV]", value=DEFAULT_ENERGY_HIGH,
#                        start=10, end=100, step=1)
#energy_step = Slider(title="energy resolution [keV]",
#                        value=DEFAULT_ENERGY_RESOLUTION, start=0.1, end=2,
#                        step=0.1)

# materials in the path
material_input = TextInput(title="Material", value=DEFAULT_MATERIAL[0])
thickness_input = TextInput(title="thickness", value=str(DEFAULT_THICKNESS[0]))
density_input = TextInput(title="density", value=str(this_material.density))

air_thickness_input = TextInput(title="air path length", value=str(DEFAULT_AIR_THICKNESS))
air_pressure_input = TextInput(title='air pressure', value=str(DEFAULT_AIR_PRESSURE))
air_temperature_input = TextInput(title='air temperature', value=str(DEFAULT_AIR_TEMPERATURE))

detector_material_input = TextInput(title='Detector', value=DEFAULT_DETECTOR_MATERIAL)
detector_thickness_input = TextInput(title='thickness', value=str(DEFAULT_DETECTOR_THICKNESS))
detector_density_input = TextInput(title='density', value=str(DEFAULT_DETECTOR_DENSITY))

p = Paragraph(text="Nothing to see here", width=500)
p.text = "All available materials: " + ", ".join(all_materials)

columns = [
    TableColumn(field="x", title="energy [keV]"),
    TableColumn(field="y", title="Percent"),
]
data_table = DataTable(source=source, columns=columns, width=400, height=700)

# the download button
download_button = Button(label="Download", button_type="success")
download_button.callback = CustomJS(args=dict(source=source),
                                    code=open(join(dirname(__file__), "download.js")).read())

log_axis_enabled = False


def update_data(attrname, old, new):

    energy = u.Quantity(
        np.arange(float(energy_low.value), float(energy_high.value),
                  float(energy_step.value)), "keV")
    x = energy.value
    i = 0

    material_list = material_input.value.split(',')
    plot_title = ""

    y = np.ones_like(x)

    if not material_input.disabled:
        for this_material_str in material_list:
            this_material = Material(str(this_material_str),
                                     u.Quantity(thickness_input.value).to("micron"),
                                     density=u.Quantity(density_input.value).to("g / cm ** 3"))
            if i == 0:
                y *= this_material.transmission(energy)
            else:
                y *= this_material.transmission(energy)
            i += 0
            plot_title += f"{this_material.name} {this_material.thickness}"

    if not air_pressure_input.disabled:
        # parse atm input
        if air_pressure_input.value.count('atm') > 0:
            air_pressure = air_pressure_input.value.split('atm')[0] * const.atm
        else:
            air_pressure = u.Quantity(air_pressure_input.value)
        air_density = get_air_density(air_pressure,
                                      u.Quantity(air_temperature_input.value))
        air = Material('air', u.Quantity(air_thickness_input.value).to("meter"),
                       density=air_density)
        y *= air.transmission(energy)
        plot_title += f" {air.name} {air.density.to('g / cm**3'):.2E} {air.thickness:.2f}"

    if not detector_material_input.disabled:
        this_detector = Material(detector_material_input.value, u.Quantity(detector_thickness_input.value),
                                 density=u.Quantity(detector_density_input.value))
        y *= this_detector.absorption(energy)
        plot_title += f" {this_detector.name} {this_detector.thickness:.2f}"

    if attrname is "log" and new:
        y = np.log10(y)
        plot.y_range.start = -4
        plot.y_range.end = 0
        plot.yaxis.axis_label = 'log(fraction)'
    else:
        plot.y_range.start = 0
        plot.y_range.end = 1
        plot.yaxis.axis_label = 'fraction'

    plot.x_range.start = energy_low.value
    plot.x_range.end = energy_high.value

    plot.title.text = plot_title
    source.data = dict(x=x, y=y)


def update_detector_density_and_data(attrname, old, new):
    detector_density_input.value = str(get_density(detector_material_input.value))
    update_data(attrname, old, new)
    return attrname, old, new

# plot.x_range.on_change('start', update_data)
# plot.x_range.on_change('end', update_data)

update_input_list = [material_input, thickness_input, density_input, air_thickness_input,
                     air_pressure_input, air_temperature_input, detector_density_input,
                     detector_thickness_input, detector_material_input,
                     #energy_low, energy_step, energy_high
                     ]

for w in update_input_list:
    w.on_change("value", update_data)

detector_material_input.on_change("value", update_detector_density_and_data)
#material_input.on_change("value", update_material_fields)

density_input.disabled = True
detector_density_input.disabled = True


def toggle_active_air(new):
    air_pressure_input.disabled = not air_pressure_input.disabled
    air_thickness_input.disabled = not air_density.disabled
    air_temperature_input.disabled = not air_temperature_input.disabled
    return new


def toggle_active(new):
    if 0 in new:
        material_input.disabled = False
        thickness_input.disabled = False
    if 0 not in new:
        material_input.disabled = True
        thickness_input.disabled = True
    if 1 in new:
        air_pressure_input.disabled = False
        air_thickness_input.disabled = False
        air_temperature_input.disabled = False
    if 1 not in new:
        air_pressure_input.disabled = True
        air_thickness_input.disabled = True
        air_temperature_input.disabled = True
    if 2 in new:
        detector_material_input.disabled = False
        detector_thickness_input.disabled = False
        detector_density_input.disabled = False
    if 2 not in new:
        detector_material_input.disabled = True
        detector_thickness_input.disabled = True
        detector_density_input.disabled = True
    update_data("toggle", 0, 0)
    return new

checkbox_group = CheckboxGroup(labels=["Material", "Air", "Detector"],
                               active=[0, 1, 2])
checkbox_group.on_click(toggle_active)

def update_plot():
    update_data("update", 0, 0)

update_plot_button = Button(label="Update", button_type="success")
update_plot_button.on_click(update_plot)


def toggle_log(new):
    if 0 in new:
        log_axis_enabled = True
    else:
        log_axis_enabled = False
    update_data("log", 0, log_axis_enabled)

plot_checkbox_group = CheckboxGroup(labels=["ylog"], active=[])
plot_checkbox_group.on_click(toggle_log)

curdoc().add_root(
    layout(
        [
            [checkbox_group, plot_checkbox_group],
            [material_input, thickness_input, density_input],
            [air_pressure_input, air_thickness_input, air_temperature_input],
            [detector_material_input, detector_thickness_input, detector_density_input],
            #[energy_low, energy_high, energy_step],
            [update_plot_button],
            [plot],
            [download_button, p]
        ],
        sizing_mode="fixed",
    )
)
curdoc().title = "Roentgen"
