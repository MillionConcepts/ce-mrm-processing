from pathlib import Path
import re
from types import MappingProxyType as MPt
import warnings

import numpy as np
import pdr
from hostess.directory import index_breadth_first
from hostess.monitors import Stopwatch

from mrm.converter.basics import modification_date

from mrm.converter.converter import PDSVersionConverter
from mrm.converter.fits import FitsLabelWriter
from mrm.shared.console import print_inline
from mrm.shared.constants import CHANNEL_FREQ, LUNAR_RADIUS_KM

warnings.simplefilter("ignore", pdr.pdr.DuplicateKeyWarning)

DISP_BOILERPLATE = """<disp:Display_Settings>
 {li_reference}
 <disp:Display_Direction>
  <disp:horizontal_display_axis>Sample</disp:horizontal_display_axis>
  <disp:horizontal_display_direction>Left to Right</disp:horizontal_display_direction>
  <disp:vertical_display_axis>Line</disp:vertical_display_axis>
  <disp:vertical_display_direction>Top to Bottom</disp:vertical_display_direction>
 </disp:Display_Direction>
</disp:Display_Settings>"""


def li_reference(identifier, reference_type):
    return f"""<Local_Internal_Reference>
  <local_identifier_reference>{identifier}</local_identifier_reference>
  <local_reference_type>{reference_type}</local_reference_type>
</Local_Internal_Reference>"""


def render_disp_boilerplate(identifier):
    return DISP_BOILERPLATE.replace(
        "{li_reference}", li_reference(identifier, "display_settings_to_array")
    )


EXTLTST_PAT = re.compile(r"(\d{1,2})_(\d{1,2})$")
BROWSE_FILEPAT = re.compile(
    r"(?P<orbiter>ce\d)?_?(?P<channel>t\d)_(?P<ptype>\w+)_"
    r"(?P<latbin>\d+_\d+)(?P<cmap>\w+)?\.png"
)


def extname_to_ltst_bin_name(extname):
    try:
        return "-".join(
            [f"{x.zfill(2)}:00" for x in EXTLTST_PAT.search(extname).groups()]
        )
    except AttributeError:
        return None


DECONV_PARAMETER_DICT = {
    "orbiter": "ORBITER",
    "channel": "CHANNEL",
    "ppd": "PPD",
    "mpp": ("PPD", lambda x: LUNAR_RADIUS_KM * 2 * np.pi / 360 / x * 1000),
}
TEMPLATE_PATH = Path(__file__).parent.parent / "label_templates"
STUB_PATH = TEMPLATE_PATH / "stubs"
BROWSE_PATH = Path(__file__).parent.parent / "bundle/browse"
MAP_PATH = Path(__file__).parent.parent / "bundle/data/maps"
LOGICAL_COORDINATE_AXES = {"LATITUDE": "y", "LONGITUDE": "x"}
PHYSICAL_COORDINATE_AXES = {"LATITUDE": "Line", "LONGITUDE": "Sample"}
INVERSE_LOGICAL_AXES = {"LATITUDE": "x", "LONGITUDE": "y"}
# TODO, maybe: we're just throwing this at all the HDUs even though not all
#  parameters are relevant for all of them. But :shrug:? Also all this dict-
#  getting is messy -- we sort of want some kind of unpacking operator, but
#  I specifically didn't want a full-on template language...
HDU_PARAMETER_DICT = {
    "ltst_bin_name": ("EXTNAME", extname_to_ltst_bin_name),
    "coordinate": ("EXTNAME", str.capitalize),
    "logical_axis": ("EXTNAME", lambda n: LOGICAL_COORDINATE_AXES.get(n)),
    "physical_axis": ("EXTNAME", lambda n: PHYSICAL_COORDINATE_AXES.get(n)),
    "other_logical_axis": ("EXTNAME", lambda n: INVERSE_LOGICAL_AXES.get(n)),
    "temp_hdu_name": (
        "EXTNAME", lambda n: re.sub(r"STDEV|LATSHIFT", "TEMP", n)
    ),
}
MAP_HDUS = ("LATSHIFT", "TEMP", "STDEV", "TBMOD", "DATMINUS")
COORD_HDUS = ("LATITUDE", "LONGITUDE")
MAPTYPES = ("temp", "latshift", "datminus", "tbmod")
# redundant, but, you know, generality
HDU_PARAMETER_DICTS = {k: HDU_PARAMETER_DICT for k in MAP_HDUS + COORD_HDUS}


MAPTYPE_NAMES = {
    "temp": "Brightness Temperature",
    "latshift": "Latitudinally-Adjusted Brightness Temperature",
    "datminus": "Model-Subtracted Brightness Temperature",
    "tbmod": "Modeled Brightness Temperature",
}


def calculate_upperleft(mpp, lines, samples):
    l0 = 1/2 - lines/2
    s0 = samples/2 - 1/2
    upleft_x = -s0 * mpp
    upleft_y = -l0 * mpp
    return upleft_x, upleft_y


def hdu_description_stub_paths():
    stubs = {p: f"{p.lower()}_hdu_description" for p in MAP_HDUS}
    stubs |= {p: f"coord_hdu_description" for p in COORD_HDUS}
    return {k: STUB_PATH / v for k, v in stubs.items()}


def map_description_stub_paths():
    return {t: STUB_PATH / f"{t}_map_description" for t in MAPTYPES}


class MapLabelWriter(FitsLabelWriter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.maptype = self.data.metaget_("PTYPE", "").lower()
        self.maptype = "temp" if self.maptype == "stdev" else self.maptype
        try:
            self.description_stub = map_description_stub_paths()[self.maptype]
        except KeyError:
            raise ValueError(f"{self.data.filename} product type unknown.")
        if self.maptype == "tbmod":
            self.template = TEMPLATE_PATH / "tbmod_template.xml"
        else:
            self.template = TEMPLATE_PATH / "map_template.xml"

    def _render_discipline_boilerplate(self):
        lrefs = {'disp_settings': [], 'cart_refs': []}
        for ident in filter(
            lambda k: k.split("_")[0] in MAP_HDUS, self.data.index
        ):
            lrefs['disp_settings'].append(render_disp_boilerplate(ident))
            lrefs['cart_refs'].append(
                li_reference(ident, "cartography_parameters_to_image_object")
            )
        return {k: "\n".join(v) for k, v in lrefs.items()}

    def _make_associations(self, reset=False, join=False):
        super()._make_associations(reset, join)
        get = self.data.metaget_
        if (timebins := re.search(r"(\d)_(\d)", get("EXTNAME"))) is None:
            raise ValueError("couldn't parse formatted time bins from name.")
        stop, start = map(int, timebins.groups())
        self.associations["ltst_binsize"] = start - stop
        self.associations |= self._render_discipline_boilerplate()
        if self.maptype != "tbmod":
            self.associations["orbiter_num"] = get("ORBITER")[-1]
        self.associations["channel"] = get("CHANNEL")
        self.associations["channel_ghz"] = CHANNEL_FREQ[get("CHANNEL")]
        self.associations["maptype"] = self.maptype
        # NOTE: assumes symmetrical +/- lat bins
        self.associations["maxlat"] = self.data["LATITUDE"].max()
        self.deletion_targets.append(self.maptype)
        self.associations["maptype_name"] = MAPTYPE_NAMES[self.maptype]
        self.associations['description'] = self._render_description()
        mapasc = [
            v for k, v in self.hdu_cvts.items()
            if k not in ("PRIMARY", "LATITUDE", "LONGITUDE")
        ][0].associations
        upx, upy = calculate_upperleft(
            mapasc["mpp"], mapasc["lines"], mapasc["samples"]
        )
        self.associations["ul_pix_x"] = upx
        self.associations["ul_pix_y"] = upy

    parameter_dicts = (DECONV_PARAMETER_DICT,)
    hdu_description_stubs = hdu_description_stub_paths()
    hdu_parameter_dicts = HDU_PARAMETER_DICTS
    description_stub: Path


# TODO, maybe: I have hardcoded pixel resolution in the text descriptions for
#  convenience. In the future we may wish to make that dynamic. If we ever
#  want to have multiple resolutions, this needs to be
#  flexible and clarified in some non-ugly way. I don't like putting the pixel
#  resolution in the browse product filenames because they're downsampled for
#  display, which makes an up-front matching ppd number misleading.


# noinspection PyMissingConstructor
class MapBrowseLabelWriter(PDSVersionConverter):
    def __init__(self, path):
        self.path = path
        self.data = None

    def _make_associations(self, reset=False, join=True):
        super()._make_associations()
        self.associations |= BROWSE_FILEPAT.search(self.path.name).groupdict()
        if self.associations['cmap'] is None:
            self.associations['cmap'] = ''
        stub = (
            f"{self.associations['ptype']}{self.associations['cmap']}"
            f"_cmap_description"
        ).lower()
        self.associations['hdu_name'] = (
            f"{self.associations['ptype']}_{self.associations['latbin']}"
        ).upper()
        self.associations['cmap_desc'] = (STUB_PATH / stub).open().read()
        self.associations['filename'] = self.path.name
        self.associations['lid_suffix'] = self.path.stem.lower()
        self.output_stem = self.path.name.replace(".png", "")
        if self.associations['ptype'] == 'latshift':
            self.associations['maptype'] = 'latshift'
        else:
            self.associations['maptype'] = 'temp'
        self.associations['modification_date'] = modification_date()
        if self.associations['ptype'] == 'tbmod':
            self.deletion_targets = ['tbmod']
        else:
            self.deletion_targets = ['change']

    template = TEMPLATE_PATH / "map_browse_template.xml"
    parameter_dicts = (MPt({}),)


if __name__ == "__main__":
    watch = Stopwatch()
    watch.start()
    for mapfile in filter(lambda p: p.suffix == ".fits", MAP_PATH.iterdir()):
        print_inline(f"making label for {mapfile.name}...")
        datawriter = MapLabelWriter(mapfile)
        datawriter.convert_label()
        datawriter.write_label(MAP_PATH)
    browse_files = index_breadth_first(BROWSE_PATH)
    for file in filter(lambda f: f['path'].endswith('.png'), browse_files):
        path = Path(file['path'])
        print_inline(f"making label for {path.name}...")
        browsewriter = MapBrowseLabelWriter(path)
        browsewriter.convert_label()
        browsewriter.write_label(path.parent)
    print(f"\n...done. ({watch.clickpeek()})")
