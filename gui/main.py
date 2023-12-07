"""
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show gui
in your browser.
"""
from pathlib import Path
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import layout, Spacer
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models import HoverTool
from bokeh.models.widgets.inputs import NumericInput
from bokeh.models.widgets import (
    Select,
    Paragraph,
    AutocompleteInput,
    CheckboxGroup,
    TableColumn,
    DataTable,
    Button,
    CheckboxGroup,
    Select,
)
from bokeh.plotting import figure
from bokeh.events import ButtonClick

import astropy.units as u
from astropy import constants as const
from astropy.units.imperial import deg_F, inch, foot, mil

from roentgen.absorption import Material, Response
from roentgen.util import get_material_density, density_ideal_gas
import roentgen

u.imperial.enable()

DEFAULT_MATERIAL = "beryllium"
DEFAULT_THICKNESS_UM = 500.0

DEFAULT_ENERGY_LOW = 1.0
DEFAULT_ENERGY_HIGH = 50.0
DEFAULT_ENERGY_STEP = 0.1

DEFAULT_DETECTOR_MATERIAL = "silicon"
DEFAULT_DETECTOR_THICKNESS_MM = 0.5
DEFAULT_DETECTOR_DENSITY = get_material_density(DEFAULT_DETECTOR_MATERIAL).value

DEFAULT_AIR_THICKNESS_M = 0.1
DEFAULT_AIR_PRESSURE = 1
DEFAULT_AIR_TEMPERATURE = 20

PLOT_HEIGHT = 400
PLOT_WIDTH = 900
TOOLS = "pan,wheel_zoom,box_zoom,box_select,undo,redo,save,reset"

custom_hover = HoverTool(
    tooltips=[
        ("energy [keV]", "@{x}{0.2f}"),
        ("transmission", "@{y}{0.4f}"),  # use @{ } for field names with spaces
    ],
    # display a tooltip whenever the cursor is vertically in line with a glyph
    mode="vline",
)

this_material = Material(DEFAULT_MATERIAL, DEFAULT_THICKNESS_UM * u.micron)
air_density = density_ideal_gas(
    DEFAULT_AIR_PRESSURE * const.atm, DEFAULT_AIR_TEMPERATURE * u.Celsius
)
air = Material("air", DEFAULT_AIR_THICKNESS_M * u.m, density=air_density)
this_detector = Material(DEFAULT_DETECTOR_MATERIAL, DEFAULT_DETECTOR_THICKNESS_MM * u.mm)

response = Response(optical_path=this_material + air, detector=this_detector)

energy = u.Quantity(np.arange(DEFAULT_ENERGY_LOW, DEFAULT_ENERGY_HIGH, DEFAULT_ENERGY_STEP), "keV")

x = energy.value
y = response.response(energy)
source = ColumnDataSource(data={"x": x, "y": y})

all_materials = list(roentgen.elements["name"]) + list(roentgen.compounds["symbol"])
all_materials.sort()
all_materials = [this_material.lower() for this_material in all_materials]

# Set up the plot
plot = figure(
    height=PLOT_HEIGHT,
    width=PLOT_WIDTH,
    tools=TOOLS,
    x_range=[1, 50],
    y_range=[0.001, 1],
)
plot.yaxis.axis_label = "Transmission fraction"
plot.xaxis.axis_label = "Energy [keV]"
plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)
plot.add_tools(custom_hover)
# Set up the inputs
ylog_checkbox = CheckboxGroup(labels=["y-log"], active=[0])

# materials in the path
material_input = AutocompleteInput(title="Material (lowercase)", value=DEFAULT_MATERIAL)
material_input.completions = all_materials
material_thickness_input = NumericInput(title="thickness", value=DEFAULT_THICKNESS_UM, mode="float")
material_density_input = NumericInput(
    title="density", value=this_material.density.value, mode="float"
)

air_thickness_input = NumericInput(
    title="air path length", value=DEFAULT_AIR_THICKNESS_M, mode="float"
)
air_pressure_input = NumericInput(title="air pressure", value=DEFAULT_AIR_PRESSURE, mode="float")
air_temperature_input = NumericInput(
    title="air temperature", value=DEFAULT_AIR_TEMPERATURE, mode="float"
)

detector_material_input = AutocompleteInput(title="Detector", value=DEFAULT_DETECTOR_MATERIAL)
detector_material_input.completions = all_materials
detector_thickness_input = NumericInput(
    title="thickness", value=DEFAULT_DETECTOR_THICKNESS_MM, mode="float"
)
detector_density_input = NumericInput(title="density", value=DEFAULT_DETECTOR_DENSITY, mode="float")

energy_high_input = NumericInput(
    title="High energy in keV", value=DEFAULT_ENERGY_HIGH, mode="float", high=2e6
)
energy_low_input = NumericInput(
    title="Low Energy in keV", value=DEFAULT_ENERGY_LOW, mode="float", low=1.0
)
energy_step_input = NumericInput(
    title="Energy step in keV", value=DEFAULT_ENERGY_STEP, mode="float"
)

p = Paragraph(text="", width=500)
p.text = f"Running roentgen version {roentgen.__version__}"

plot_title = Paragraph(text="", width=500)
plot_title.text = f"{response}"

columns = [
    TableColumn(field="x", title="energy [keV]"),
    TableColumn(field="y", title="Percent"),
]
data_table = DataTable(source=source, columns=columns, width=400, height=700)

download_button = Button(label="Download", button_type="success")
download_button.js_on_event(
    "button_click",
    CustomJS(
        args=dict(source=source),
        code=(Path(__file__).parent / "download.js").read_text("utf8"),
    ),
)


def convert_air_pressure(value, current_unit, new_unit):
    if current_unit == "atm":
        air_pressure = u.Quantity(value * const.atm, "Pa")
    elif current_unit == "torr":
        air_pressure = u.Quantity(value * const.atm / 760.0, "Pa")
    else:
        air_pressure = u.Quantity(value, current_unit)

    if new_unit == "atm":
        return (air_pressure.to("Pa") / const.atm).value
    elif new_unit == "torr":
        return (air_pressure.to("Pa") / const.atm).value * 760.0
    else:
        return air_pressure.to(new_unit)


def update_response(attrname, old, new):
    # check whether the input variables have changed and update the response
    global response

    if not material_input.disabled:
        if material_input.value.lower() in all_materials:
            this_thickness = u.Quantity(material_thickness_input.value, material_thick_unit.value)
            this_density = u.Quantity(material_density_input.value, material_density_unit.value)
            this_material = Material(
                material_input.value.lower(), this_thickness, density=this_density
            )
    else:
        # if material not selected, just make a bogus material with no thickness
        this_material = Material(material_input.value.lower(), 0 * u.mm)

    if not air_pressure_input.disabled:
        if air_pressure_unit.value == "atm":
            air_pressure = u.Quantity(air_pressure_input.value * const.atm, "Pa")
        elif air_pressure_unit.value == "torr":
            air_pressure = u.Quantity(air_pressure_input.value * const.atm / 760.0, "Pa")
        else:
            air_pressure = u.Quantity(air_pressure_input.value, air_pressure_unit.value)

        air_path_length = u.Quantity(air_thickness_input.value, air_thick_unit.value)
        air_temperature = u.Quantity(air_temperature_input.value, air_temp_unit.value).to(
            "Celsius", equivalencies=u.temperature()
        )
        air_density = density_ideal_gas(air_pressure, air_temperature)
        air = Material("air", air_path_length, density=air_density)
    else:
        # if air is not selected than just add bogus air with no thickness
        air = Material("air", 0 * u.mm, density=0 * u.g / u.cm**3)

    if not detector_material_input.disabled:
        if detector_material_input.value.lower() in all_materials:
            this_thickness = u.Quantity(detector_thickness_input.value, detector_thick_unit.value)
            this_density = u.Quantity(detector_density_input.value, detector_density_unit.value)
            this_detector = Material(
                detector_material_input.value.lower(), this_thickness, density=this_density
            )
    else:
        # if detector is not selected than add bogus super thick detector
        this_detector = Material(detector_material_input.value.lower(), 1 * u.Mm)
    response = Response(optical_path=this_material + air, detector=this_detector)


def update_data(attrname, old, new):
    global source
    energy = u.Quantity(
        np.arange(energy_low_input.value, energy_high_input.value, energy_step_input.value), "keV"
    )
    x = energy.value
    y = response.response(energy)

    plot_title.text = f"{response}"

    if plot_checkbox_group.active:
        y = np.log10(response.response(energy))
        plot.y_range.start = -4
        plot.y_range.end = np.max(y)
        if not detector_material_input.disabled:
            plot.yaxis.axis_label = "log(Response fraction)"
        else:
            plot.yaxis.axis_label = "log(Transmission fraction)"
    else:
        plot.y_range.start = 0
        plot.y_range.end = 1
        if not detector_material_input.disabled:
            plot.yaxis.axis_label = "Response fraction)"
        else:
            plot.yaxis.axis_label = "Transmission fraction"

    source.data = dict(x=x, y=y)


def toggle_active(attr, old, new):
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


checkbox_group = CheckboxGroup(labels=["Material", "Air", "Detector"], active=[0, 1, 2])
checkbox_group.on_change("active", toggle_active)


def update_button_action():
    update_response("update_plot_button", 0, 0)
    update_data("update", 0, 0)


update_plot_button = Button(label="Update Plot", button_type="success")
update_plot_button.on_click(update_button_action)

# plot.x_range.on_change("start", update_data)
# plot.x_range.on_change("end", update_data)


plot_checkbox_group = CheckboxGroup(labels=["ylog"], active=[])

length_units = ["m", "mm", "micron", "inch", "foot", "mil"]
pressure_units = ["Pa", "torr", "atm"]
density_units = ["g / cm ** 3", "kg / m ** 3"]
temperature_units = ["K", "deg_C", "deg_F"]


material_thick_unit = Select(title="unit", value=length_units[2], options=length_units)


def update_mat_thick_units(attr, old, new):
    material_thickness_input.value = u.Quantity(material_thickness_input.value, old).to(new).value


material_thick_unit.on_change("value", update_mat_thick_units)


detector_thick_unit = Select(title="unit", value=length_units[1], options=length_units)


def update_det_thick_units(attr, old, new):
    detector_thickness_input.value = u.Quantity(detector_thickness_input.value, old).to(new).value


detector_thick_unit.on_change("value", update_det_thick_units)
air_thick_unit = Select(title="unit", value=length_units[0], options=length_units)


def update_air_thick_units(attr, old, new):
    air_thickness_input.value = u.Quantity(air_thickness_input.value, old).to(new).value


air_thick_unit.on_change("value", update_air_thick_units)

material_density_unit = Select(title="unit", value=density_units[0], options=density_units)


def update_mat_density_units(attr, old, new):
    material_density_input.value = u.Quantity(material_density_input.value, old).to(new).value


material_density_unit.on_change("value", update_mat_density_units)

detector_density_unit = Select(title="unit", value=density_units[0], options=density_units)


def update_det_density_units(attr, old, new):
    detector_density_input.value = u.Quantity(detector_density_input.value, old).to(new).value


detector_density_unit.on_change("value", update_det_density_units)
air_pressure_unit = Select(title="unit", value=pressure_units[2], options=pressure_units)


def update_air_pressure_units(attr, old, new):
    air_pressure = convert_air_pressure(air_pressure_input.value, old, new)
    air_pressure_input.value = air_pressure


air_pressure_unit.on_change("value", update_air_pressure_units)

air_temp_unit = Select(title="unit", value=temperature_units[1], options=temperature_units)


def update_air_temp_units(attr, old, new):
    air_temperature_input.value = (
        u.Quantity(air_temperature_input.value, old).to(new, equivalencies=u.temperature()).value
    )


air_temp_unit.on_change("value", update_air_temp_units)


def update_material_density(attr, old, new):
    # if the material is changed then update the density
    material_density_input.value = get_material_density(material_input.value).value


material_input.on_change("value", update_material_density)


def update_detector_density(attr, old, new):
    # if the material is changed then update the density
    detector_density_input.value = get_material_density(detector_material_input.value).value


detector_material_input.on_change("value", update_detector_density)

curdoc().add_root(
    layout(
        [
            [plot_title],
            [plot],
            [checkbox_group, plot_checkbox_group],
            [energy_low_input, energy_high_input, energy_step_input],
            [
                material_input,
                material_thickness_input,
                material_thick_unit,
                material_density_input,
                material_density_unit,
            ],
            [
                air_pressure_input,
                air_pressure_unit,
                air_thickness_input,
                air_thick_unit,
                air_temperature_input,
                air_temp_unit,
            ],
            [
                detector_material_input,
                detector_thickness_input,
                detector_thick_unit,
                detector_density_input,
                detector_density_unit,
            ],
            # [energy_low, energy_high, energy_step],
            [download_button, Spacer(), update_plot_button],
            [p],
        ],
        sizing_mode="scale_width",
    )
)
curdoc().title = "Roentgen"
