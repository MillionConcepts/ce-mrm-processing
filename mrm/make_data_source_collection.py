from itertools import chain
from pathlib import Path
import re
import shutil

from dateutil.parser import parse as dtp
from dustgoggles.structures import MaybePool

from mrm.converter.basics import modification_date
from mrm.converter.converter import PDSVersionConverter

L2C_PARAMETER_DICT = {
    'orbiter_num': ('SPACECRAFT_NAME', lambda n: n[-1]),
    'product_creation_time': 'PRODUCT_CREATION_TIME',
    'rows': 'ROWS'
}
ORBITPAT = re.compile(r".*_(\d{4})_[AB]\.2[BC]")

TEMPLATE_PATH = Path(__file__).parent.parent / "label_templates"
STUB_PATH = TEMPLATE_PATH / "stubs"
OUTPUT_ROOT = Path(__file__).parent.parent / "bundle" / "data_source"
CE1_L2C_PATH = Path(
    "/datascratch/CE1_MRM_2C_20071127161816_20090114192053_B/DATA"
)
CE2_L2C_PATH = Path(
    "/datascratch/CE2_MRM_2C_20101015085002_20110520125150_A/DATA"
)


class L2CLabelMaker(PDSVersionConverter):

    def make_associations(self, reset=True, join=True):
        super()._make_associations(reset, join)
        get, assoc = self.data.metaget_, self.associations
        num = assoc['orbiter_num']
        assoc['filename'] = Path(self.data.filename).name
        assoc['filename_stem'] = Path(self.data.filename).stem
        assoc['lid_suffix'] = Path(self.data.filename).stem.lower()
        with (STUB_PATH / f'ce{num}_l2c_external_reference').open() as stream:
            assoc['external_reference'] = stream.read()
        assoc['header_size'] = get('RECORD_BYTES') * get('LABEL_RECORDS')
        with (STUB_PATH / f'ce{num}_l2c_record_character').open() as stream:
            assoc['record_character'] = stream.read()
        assoc['orbit_number'] = int(
            ORBITPAT.search(self.data.filename).group(1)
        )
        assoc['start_date_time'] = self.data.TABLE['Time'][0]
        assoc['stop_date_time'] = self.data.TABLE['Time'].iloc[-1]
        assoc['modification_date'] = modification_date()

    template = TEMPLATE_PATH / "l2c_template.xml"
    parameter_dicts = (L2C_PARAMETER_DICT,)


def subdir(llm: L2CLabelMaker):
    orbiter = f"ce{llm.associations['orbiter_num']}"
    sdtime = dtp(llm.associations['start_date_time'])
    return Path(orbiter) / f"{sdtime.year}{str(sdtime.month).zfill(2)}"


def make_l2c_product(source_path):
    cvt = L2CLabelMaker(source_path)
    # we're not including 'empty' tables
    if cvt.data.metaget_("ROWS") == 0:
        return
    cvt.make_associations()
    outpath = OUTPUT_ROOT / subdir(cvt)
    outpath.mkdir(parents=True, exist_ok=True)
    cvt.convert_label()
    cvt.write_label(outpath)
    shutil.copy(source_path, outpath)


def make_l2c_products(source_paths):
    for sp in source_paths:
        make_l2c_product(sp)


THREADS = 8

if __name__ == "__main__":
    shutil.rmtree(OUTPUT_ROOT / "ce1", ignore_errors=True)
    shutil.rmtree(OUTPUT_ROOT / "ce2", ignore_errors=True)
    argrecs = [
        {'args': (sp,)}
        for sp in chain(*(CE1_L2C_PATH.iterdir(), CE2_L2C_PATH.iterdir()))
    ]
    pool = MaybePool(THREADS)
    pool.map(make_l2c_product, argrecs, as_chunks=True)
    pool.close()
    pool.join()
