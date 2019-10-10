"""
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show gui
in your browser.
"""
import numpy as np
from os.path import dirname, join

from bokeh.io import curdoc
from bokeh.layouts import layout, Spacer
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
    Select,
    Div
)
from bokeh.plotting import figure

import astropy.units as u
from astropy import constants as const
from astropy.units.imperial import deg_F, inch, foot, mil
u.imperial.enable()

from roentgen.absorption import Material, get_density
import roentgen

R = 287.058 * u.J / u.kg / u.Kelvin


def get_air_density(pressure, temperature):
    return pressure / (R * temperature.to('K', equivalencies=u.temperature()))

DEFAULT_MATERIAL = ["silicon"]
DEFAULT_THICKNESS = [100.0]
DEFAULT_ENERGY_LOW = 1.0
DEFAULT_ENERGY_HIGH = 200.0
DEFAULT_DENSITY = [2.33]
NUMBER_OF_MATERIALS = len(DEFAULT_MATERIAL)

DEFAULT_DETECTOR_MATERIAL = 'cdte'
DEFAULT_DETECTOR_THICKNESS = 1
DEFAULT_DETECTOR_DENSITY = get_density(DEFAULT_DETECTOR_MATERIAL).value

DEFAULT_ENERGY_RESOLUTION = 0.25

DEFAULT_AIR_THICKNESS = 1
DEFAULT_AIR_PRESSURE = 1
DEFAULT_AIR_TEMPERATURE = 25

PLOT_HEIGHT = 300
PLOT_WIDTH = 900
TOOLS = "pan,box_zoom,box_select,undo,redo,save,reset"

custom_hover = HoverTool(
    tooltips=[
        ('energy [keV]',   '@{x}{0.2f}'            ),
        ('transmission',  '@{y}{0.2f}' ), # use @{ } for field names with spaces
    ],

    # display a tooltip whenever the cursor is vertically in line with a glyph
    mode='vline'
)

# defaults
material_list = []
energy = u.Quantity(np.arange(DEFAULT_ENERGY_LOW, DEFAULT_ENERGY_HIGH,
                              DEFAULT_ENERGY_RESOLUTION), "keV")

this_material = Material(DEFAULT_MATERIAL[0], DEFAULT_THICKNESS[0] * u.micron)
y = this_material.transmission(energy)
air_density = get_air_density(DEFAULT_AIR_PRESSURE * const.atm, DEFAULT_AIR_TEMPERATURE * u.Celsius)
air = Material('air', DEFAULT_AIR_THICKNESS * u.m, density=air_density)
y *= air.transmission(energy)
this_detector = Material(DEFAULT_DETECTOR_MATERIAL, DEFAULT_DETECTOR_THICKNESS * u.mm)
y *= this_detector.absorption(energy)

x = energy.value

source = ColumnDataSource(data={"x": x, "y": y})

all_materials = list(roentgen.elements["name"]) + \
                list(roentgen.compounds["symbol"])
all_materials.sort()
all_materials = [this_material.lower() for this_material in all_materials]

# Set up the plot
plot = figure(
    plot_height=PLOT_HEIGHT,
    plot_width=PLOT_WIDTH,
    tools=TOOLS,
    x_range=[1, 50],
    y_range=[0, 1],
)
plot.yaxis.axis_label = "Transmission fraction"
plot.xaxis.axis_label = "Energy [keV]"
plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)
plot.title.text = "Plot Title"
plot.add_tools(custom_hover)
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
material_input = AutocompleteInput(title="Material (lowercase)", value=DEFAULT_MATERIAL[0])
material_input.completions = all_materials
material_thickness_input = TextInput(title="thickness", value=str(DEFAULT_THICKNESS[0]))
material_density_input = TextInput(title="density", value=str(this_material.density.value))

air_thickness_input = TextInput(title="air path length", value=str(DEFAULT_AIR_THICKNESS))
air_pressure_input = TextInput(title='air pressure', value=str(DEFAULT_AIR_PRESSURE))
air_temperature_input = TextInput(title='air temperature', value=str(DEFAULT_AIR_TEMPERATURE))

detector_material_input = AutocompleteInput(title='Detector', value=DEFAULT_DETECTOR_MATERIAL)
detector_material_input.completions = all_materials
detector_thickness_input = TextInput(title='thickness', value=str(DEFAULT_DETECTOR_THICKNESS))
detector_density_input = TextInput(title='density', value=str(DEFAULT_DETECTOR_DENSITY))

p = Paragraph(text="", width=500)
p.text = "Messages: "

columns = [
    TableColumn(field="x", title="energy [keV]"),
    TableColumn(field="y", title="Percent"),
]
data_table = DataTable(source=source, columns=columns, width=400, height=700)

# the download button
download_button = Button(label="Download", button_type="success")
download_button.callback = CustomJS(args=dict(source=source),
                                    code=open(join(dirname(__file__), "download.js")).read())


def convert_air_pressure(value, current_unit, new_unit):
    if current_unit == "atm":
        air_pressure = u.Quantity(value * const.atm, "Pa")
    elif current_unit == "torr":
        air_pressure = u.Quantity(value * const.atm / 760., "Pa")
    else:
        air_pressure = u.Quantity(value, current_unit)

    if new_unit == "atm":
        return (air_pressure.to("Pa") / const.atm).value
    elif new_unit == "torr":
        return (air_pressure.to("Pa") / const.atm).value * 760.0
    else:
        return air_pressure.to(new_unit)


def update_data(attrname, old, new):

    i = 0
    plot_title = ""
    p.text = 'Messages: '
    y = np.ones_like(x)

    if not material_input.disabled:
        this_material_str = material_input.value
        if str(this_material_str.lower()) in all_materials:
            this_thickness = u.Quantity(material_thickness_input.value, material_thick_unit.value)
            this_density = u.Quantity(material_density_input.value, material_density_unit.value)

            this_material = Material(this_material_str, this_thickness,
                                     density=this_density)
            if i == 0:
                y *= this_material.transmission(energy)
            else:
                y *= this_material.transmission(energy)
            i += 0
            plot_title += f"{this_material.name} {this_material.thickness}"
        else:
            p.text += f'material {this_material_str} not recognized'

    if not air_pressure_input.disabled:
        if air_pressure_unit.value == "atm":
            air_pressure = u.Quantity(air_pressure_input.value * const.atm, "Pa")
        elif air_pressure_unit.value == "torr":
            air_pressure = u.Quantity(air_pressure_input.value * const.atm / 760., "Pa")
        else:
            air_pressure = u.Quantity(air_pressure_input.value, air_pressure_unit.value)

        air_path_length = u.Quantity(air_thickness_input.value, air_thick_unit.value)
        air_temperature = u.Quantity(air_temperature_input.value,
                                     air_temp_unit.value).to("Celsius", equivalencies=u.temperature())

        # parse atm input
        #if air_pressure_input.value.count('atm') > 0:
        #    air_pressure = air_pressure_input.value.split('atm')[0] * const.atm
        #else:
        air_density = get_air_density(air_pressure, air_temperature)
        air = Material('air', air_path_length, density=air_density)
        y *= air.transmission(energy)
        plot_title += f" {air.name} {air.density.to('g / cm**3'):.2E} {air.thickness:.2f}"

    if not detector_material_input.disabled:
        this_detector_material_str = detector_material_input.value
        if this_detector_material_str.lower() in all_materials:
            this_thickness = u.Quantity(detector_thickness_input.value, detector_thick_unit.value)
            this_density = u.Quantity(detector_density_input.value, detector_density_unit.value)
            this_detector = Material(this_detector_material_str, this_thickness,
                                     density=this_density)
            y *= this_detector.absorption(energy)
            plot_title += f" {this_detector.name} {this_detector.thickness:.2f}"
        else:
            p.text += f'detector {this_detector_material_str} not recognized'

    if plot_checkbox_group.active:
        y = np.log10(y)
        plot.y_range.start = -4
        plot.y_range.end = 0
        plot.yaxis.axis_label = 'log(Transmission fraction)'
    else:
        plot.y_range.start = 0
        plot.y_range.end = 1
        plot.yaxis.axis_label = 'Transmission fraction'

    plot.title.text = plot_title
    source.data = dict(x=x, y=y)


def toggle_active(new):
    if 0 in new:
        material_input.disabled = False
        material_thick_unit.disabled = False
        material_density_input.disabled = False
    if 0 not in new:
        material_input.disabled = True
        material_thick_unit.disabled = True
        material_density_input.disabled = True
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
    return new

checkbox_group = CheckboxGroup(labels=["Material", "Air", "Detector"],
                               active=[0, 1, 2])
checkbox_group.on_click(toggle_active)


def update_plot():
    update_data("update", 0, 0)

update_plot_button = Button(label="Update Plot", button_type="success")
update_plot_button.on_click(update_plot)

plot_checkbox_group = CheckboxGroup(labels=["ylog"], active=[])

length_units = ["m", "mm", "micron", "inch", "foot", "mil"]
pressure_units = ["Pa", "torr", "atm"]
density_units = ["g / cm ** 3", "kg / m ** 3"]
temperature_units = ["K", "deg_C", "deg_F"]



material_thick_unit = Select(title="unit", value=length_units[2], options=length_units)
def update_mat_thick_units(attr, old, new):
    material_thickness_input.value = str(u.Quantity(material_thickness_input.value, old).to(new).value)
material_thick_unit.on_change('value', update_mat_thick_units)

detector_thick_unit = Select(title="unit", value=length_units[1], options=length_units)
def update_det_thick_units(attr, old, new):
    detector_thickness_input.value = str(u.Quantity(detector_thickness_input.value, old).to(new).value)
detector_thick_unit.on_change('value', update_det_thick_units)

air_thick_unit = Select(title="unit", value=length_units[0], options=length_units)
def update_air_thick_units(attr, old, new):
    air_thickness_input.value = str(u.Quantity(air_thickness_input.value, old).to(new).value)
air_thick_unit.on_change('value', update_air_thick_units)

material_density_unit = Select(title="unit", value=density_units[0], options=density_units)
def update_mat_density_units(attr, old, new):
    material_density_input.value = str(u.Quantity(material_density_input.value, old).to(new).value)
material_density_unit.on_change('value', update_mat_density_units)

detector_density_unit = Select(title="unit", value=density_units[0], options=density_units)
def update_det_density_units(attr, old, new):
    detector_density_input.value = str(u.Quantity(detector_density_input.value, old).to(new).value)
detector_density_unit.on_change('value', update_det_density_units)

air_pressure_unit = Select(title="unit", value=pressure_units[2], options=pressure_units)
def update_air_pressure_units(attr, old, new):
    air_pressure = convert_air_pressure(air_pressure_input.value, old, new)
    air_pressure_input.value = str(air_pressure)
air_pressure_unit.on_change('value', update_air_pressure_units)


air_temp_unit = Select(title="unit", value=temperature_units[1], options=temperature_units)
def update_air_temp_units(attr, old, new):
    air_temperature_input.value = str(u.Quantity(air_temperature_input.value, old).to(new, equivalencies=u.temperature()).value)
air_temp_unit.on_change('value', update_air_temp_units)


curdoc().add_root(
    layout(
        [
            [plot],
            [checkbox_group, plot_checkbox_group],
            [material_input, material_thickness_input, material_thick_unit, material_density_input, material_density_unit],
            [air_pressure_input, air_pressure_unit, air_thickness_input, air_thick_unit, air_temperature_input, air_temp_unit],
            [detector_material_input, detector_thickness_input, detector_thick_unit, detector_density_input, detector_density_unit],
            #[energy_low, energy_high, energy_step],
            [download_button, Spacer(), update_plot_button],
            [p]
        ],
        sizing_mode="scale_width",
    )
)
curdoc().title = "Roentgen"
