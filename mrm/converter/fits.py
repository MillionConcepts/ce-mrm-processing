from itertools import chain
from pathlib import Path
import re
from types import MappingProxyType as MPt, NoneType
from typing import Mapping, Optional, Union, Never

from astropy.io import fits
from cytoolz import groupby, keyfilter, valfilter
from dustgoggles.func import constant
from dustgoggles.structures import dig_and_edit
from multidict import MultiDict
import pdr
from pdr.parselabel.pds4 import log_and_pass

from mrm.converter.converter import PDSVersionConverter
from mrm.converter.templating import xmlwrap


# TODO: this would be better as a simple set of regexes
def get_by_substr[T](mapping: Mapping[str, T], key: str, default=None) -> T:
    """
    Get a value from a mapping in a slightly fuzzy fashion. Prefers exact
    match between a key of the mapping and `keyname`; if it finds none, falls
    back to case-insensitive match, then case-insensitive substring match.
    """
    if key in mapping.keys():
        return mapping[key]
    lowmap = groupby(lambda k: k.lower(), mapping.keys())
    opts = lowmap.get(key.lower())
    if opts is not None:
        return mapping[opts[0]]
    if len(
        keygroups := keyfilter(lambda k: k.lower() in key.lower(), lowmap)
    ) == 0:
        return default
    if len(keyopts := keygroups[tuple(keygroups.keys())[0]]) == 0:
        return default
    return mapping[keyopts[0]]


def imsz_from_header(header) -> tuple[int]:
    """
    get image size from either compressed or uncompressed FITS image headers.
    should work on headers returned from either fitsio or astropy,
    returns in 'reversed' order for numpy array indexing.
    """
    key_type = "ZNAXIS" if "ZNAXIS" in header.keys() else "NAXIS"
    axis_entries: dict[str, int] = keyfilter(
        lambda k: False if k is None else re.match(rf"{key_type}\d", k),
        dict(header),
    )
    # noinspection PyTypeChecker
    return tuple(reversed(axis_entries.values()))


# TODO: implement table handling
def fitsstat(obj: Union[str, Path, fits.HDUList]) -> dict:
    """
    produce a dict describing the characteristics
    of a FITS file and its extensions.
    """
    info = {}
    hdul = obj if isinstance(obj, fits.HDUList) else fits.open(obj)
    hdulinfo = hdul.info(False)
    for hdu_ix, hdu in enumerate(hdul):
        hinfo = hdu.fileinfo()
        hdu_info = {
            "hdu_size": hinfo["datLoc"] - hinfo["hdrLoc"] + hinfo["datSpan"],
            "header_size": hinfo['datLoc'] - hinfo['hdrLoc'],
            "data_size": hinfo["datSpan"],
            "header_offset": hinfo['hdrLoc'],
            "data_offset": hinfo["datLoc"],
            "name": hdulinfo[hdu_ix][1],
            "hdu_type": hdulinfo[hdu_ix][3],
            "dim": imsz_from_header(hdu.header),
        }
        if len(hdu_info["dim"]) == 0:
            hdu_info["itemsize"] = None
        else:
            hdu_info["itemsize"] = abs(hdu.header["BITPIX"])
        info[hdu_ix] = hdu_info
    return info


BITPIX_PDS4_MAPPING = MPt(
    {
        -64: "IEEE754MSBDouble",
        -32: "IEEE754MSBSingle",
        8: "UnsignedByte",
        16: "SignedMSB2",
        32: "SignedMSB4",
        64: "SignedMSB8",
    }
)

IMAGE_HDU_TAGS = MPt(
    {
        "dimensionality": "NAXIS",
        "missing_constant": "BLANK",
        "lines": "NAXIS2",
        "bands": "NAXIS3",
        "samples": "NAXIS1",
        "offset": "BZERO",
        "scale": "BSCALE",
        "dtype": ("BITPIX", lambda bitpix: BITPIX_PDS4_MAPPING[bitpix]),
        "unit": "BUNIT"
    }
)

HDU_PARAMETER_DICTS = MPt(
    {"ImageHDU": IMAGE_HDU_TAGS, "PrimaryHDU": IMAGE_HDU_TAGS}
)

STUB_TEMPLATE_PATH = Path(__file__).parent.absolute() / "fits_template_stubs"

STUB_TEMPLATES = MPt(
    {
        'ImageHDU': STUB_TEMPLATE_PATH / "image_hdu_stub.xml",
        'PrimaryHDU': STUB_TEMPLATE_PATH / "image_hdu_stub.xml",
        "header": STUB_TEMPLATE_PATH / "fits_header_stub.xml"
    }
)


def _truncate_metadata(
    metadata: pdr.Metadata, blockname: str, standard: str = "FITS"
):
    mapping, params = metadata.metablock(blockname), []
    if mapping is None:
        return MultiDict()
    dig_and_edit(
        mapping,
        constant(True),
        lambda k, v: log_and_pass(params, k, v),
        mtypes=(MultiDict,),
    )
    return pdr.Metadata((mapping, params), standard)


class HDUStubConverter(PDSVersionConverter):
    # noinspection PyMissingConstructor
    def __init__(self, *_, **__):
        raise TypeError(
            "Do not initialize this object directly. Use from_converter()."
        )

    def _render_description(self):
        stub = self.description_stub.open().read()
        description = self.fill_template(
            text=xmlwrap('description')(stub), pprint_xml=False
        )
        return re.sub(r"<\?xml.*?\?>\n", "", description)

    def _make_associations(
        self, reset: NoneType = None, join: NoneType = None
    ):
        super()._make_associations(join=True)
        if self.associations["dimensionality"] > 3:
            raise NotImplementedError(
                "FITS arrays with more than 3 dimensions are not supported."
            )
        self.associations['hdu_name'] = self.hdu_name
        self.associations["data_offset"] = self.data_offset
        self.associations["header_offset"] = self.header_offset
        self.associations["header_size"] = self.header_size
        self.deletion_targets = [
            f"no_{k}" for k in self.parameter_dicts[0].keys()
            if self.associations.get(k) is None
        ]
        self.deletion_targets.append(f"{self.associations['dimensionality']}d")
        if self.description_stub is not None:
            self.associations["description"] = self._render_description()

    def convert_label(self, pprint_xml: Never = None):
        header_xml = self.fill_template(
            STUB_TEMPLATES["header"], pprint_xml=False
        )
        self.header_xml = re.sub(r"<\?xml.*?\?>\n", "", header_xml)
        if self.data_size == 0:
            return
        data_xml = self.fill_template(
            STUB_TEMPLATES[self.hdu_type], pprint_xml=False
        )
        self.data_xml = re.sub(r"<\?xml.*?\?>\n", "", data_xml)

    def build_hdu_label(self):
        self.convert_label()
        return [self.header_xml, self.data_xml]

    @classmethod
    def from_converter(cls, hdul_cvt, hdu_name=None, **kwargs):
        if hdu_name is None:
            raise TypeError("Must pass hdu_name to create this object.")
        cvt = super().__new__(cls)
        cvt.hdu_name = hdu_name
        for k in cls.syntactic_kwargs:
            setattr(cvt, k, kwargs.get(k))
        hdu_name_matches = valfilter(
            lambda v: v["name"] == cvt.hdu_name, hdul_cvt.fitsinfo
        )
        if len(hdu_name_matches) == 0:
            raise ValueError(f"No HDU named {cvt.hdu_name}.")
        if len(hdu_name_matches) > 1:
            raise ValueError(f"Multiple HDUs named {cvt.hdu_name}.")
        cvt.hdu_ix, hdu_info = list(hdu_name_matches.items())[0]
        for k, v in filter(lambda kv: kv[0] != "name", hdu_info.items()):
            setattr(cvt, k, v)
        if cvt.hdu_type not in HDU_PARAMETER_DICTS.keys():
            raise NotImplementedError(f"{cvt.hdu_type} not yet supported.")
        cvt.parameter_dicts = [HDU_PARAMETER_DICTS[cvt.hdu_type],]
        if (extra_pdict := kwargs.pop("parameter_dict", None)) is not None:
            cvt.parameter_dicts.append(extra_pdict)
        cvt.associations = hdul_cvt.associations | {
            k: v for k, v in kwargs.items() if k not in cls.syntactic_kwargs
        }
        cvt.metadata = _truncate_metadata(hdul_cvt.metadata, cvt.hdu_name)
        cvt._make_associations()
        return cvt

    data_offset: int
    header_offset: int
    data_size: int
    description_stub: Optional[Path]
    header_size: int
    hdu_name: str
    hdu_type: str
    data_xml: Optional[str] = None
    header_xml: str
    syntactic_kwargs = ("description_stub", "unit")


class FitsLabelWriter(PDSVersionConverter):

    def __init__(self, filename: Union[str, Path]):
        self.data = pdr.fastread(filename)
        if self.data.standard != "FITS":
            raise OSError(f"{filename} does not appear to be a FITS file.")
        self.output_paths = {'fits': filename}
        super().__init__(skip_load=True)
        self.output_stem = Path(filename).stem
        self.fitsinfo = fitsstat(self.data._hdulist)
        for attr in ("hdu_description_stubs", "hdu_parameter_dicts", "units"):
            if getattr(self, attr) is None:
                setattr(self, attr, {})

    def _render_description(self):
        stub = self.description_stub.open().read()
        description = self.fill_template(
            text=xmlwrap('description')(stub), pprint_xml=False
        )
        return re.sub(r"<\?xml.*?\?>\n", "", description)

    def _make_hdu_converter(self, hdu_name):
        return HDUStubConverter.from_converter(
            self,
            hdu_name=hdu_name,
            unit=self.units.get(hdu_name),
            description_stub=get_by_substr(self.hdu_description_stubs, hdu_name),
            parameter_dict=get_by_substr(self.hdu_parameter_dicts, hdu_name)
        )

    def _make_associations(self, reset=False, join=False):
        super()._make_associations(reset, join)
        self.hdu_cvts = {
            hdu_name: self._make_hdu_converter(hdu_name)
            for hdu_name in self.data.keys()
        }
        self.hdu_labels = {
            hdu_name: hcvt.build_hdu_label()
            for hdu_name, hcvt in self.hdu_cvts.items()
        }
        # NOTE: If an element of `hdu_objects` is `None`, it means that the
        #  data section of the corresponding HDU is an array with no elements,
        #  most likely a "stub" PRIMARY HDU. We could legally label it as a
        #  0-length Array_1D, but it is cleaner to not mention it at all.
        self.associations['hdu_objects'] = "".join(
            filter(None, chain.from_iterable(self.hdu_labels.values()))
        )

    def convert_label(self, pprint_xml=True):
        self._make_associations(False, True)
        self.pds4_label = self.fill_template(pprint_xml=pprint_xml)

    description_stub = None
    hdu_parameter_dicts = None
    hdu_description_stubs = None
    units = None
