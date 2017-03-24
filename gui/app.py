import flask

import astropy.units as u
import numpy as np
import roentgen

from roentgen.absorption import Material

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
from bokeh.models import ColumnDataSource, CustomJS, HoverTool

app = flask.Flask(__name__)

DEFAULT_MATERIAL = 'Si'
DEFAULT_THICKNESS = 100.0

PLOT_HEIGHT = 300
PLOT_WIDTH = 900
TOOLS = 'pan,box_zoom,box_select,crosshair,undo,redo,save,reset'

# defaults
mat = Material(DEFAULT_MATERIAL, DEFAULT_THICKNESS * u.micron)
energy = u.Quantity(np.arange(1, 1000), 'keV')

source = ColumnDataSource(data={'x': energy.value, 'y':mat.absorption(energy).value})
source_static = ColumnDataSource(data={'x': energy.value, 'y':mat.absorption(energy).value})

@app.route('/')
def index():
    args = flask.request.args
    _input_material = str(args.get('_input_material', DEFAULT_MATERIAL))
    _input_thickness = str(args.get('_input_thickness', DEFAULT_THICKNESS))

    mat = Material(_input_material, float(_input_thickness) * u.micron)
    energy = u.Quantity(np.arange(1, 100), 'keV')

    source = ColumnDataSource(data={'x': energy.value, 'y': mat.absorption(energy).value})
    source_static = ColumnDataSource(data={'x': energy.value, 'y': mat.absorption(energy).value})

    fig = figure(title="Absorption", tools=TOOLS, plot_height=PLOT_HEIGHT, width=PLOT_WIDTH, toolbar_location="right",
                 x_axis_label='Energy [keV]')
    fig.line('x', 'y', source=source_static, color='red', line_width=2, legend=_input_material)

    hover = HoverTool()
    #hover.tooltips = [
    #    ("x", "@time_str"),
    #    ("y", "@xrsb")
    #]
    fig.add_tools(hover)

    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = components(fig, INLINE)
    html = flask.render_template(
        'index.html',
        plot_script=script,
        plot_div=div,
        material_list=list(roentgen.elements['name']),
        js_resources=js_resources,
        css_resources=css_resources,
        _input_material=_input_material,
        _input_thickness=_input_thickness,
    )
    return encode_utf8(html)

if __name__ == '__main__':
    print(__doc__)
    app.run(debug=True)
