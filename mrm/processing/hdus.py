import numpy as np
from astropy.io import fits
from dustgoggles.func import gmap


def to_scaled_hdu(arr, dtype="int16", max_abs_err=2e-3, **hdu_kwargs):
    masked = np.ma.masked_invalid(arr)
    hdu = fits.ImageHDU(masked.filled(masked.mean()), **hdu_kwargs)
    hdu.scale(dtype, "minmax")
    typemin = np.iinfo(getattr(np, dtype)).min
    if (apply_mask := bool(masked.mask.any())) is True:
        if hdu.data.min() == typemin:
            raise NotImplementedError("cannot mask; min type value used")
        hdu.data[masked.mask], hdu.header["BLANK"] = (typemin, typemin)
    del masked
    if max_abs_err is not None:
        # this little set of contortions replicates astropy's HDU-loading
        # behavior without forcing us to serialize or bounce to disk
        check = (
            (np.ma.masked_equal(hdu.data, typemin) if apply_mask else hdu.data)
            # NOTE: these keywords may not exist if the input array didn't
            # 'need' scaling / offset for whatever reason
            * hdu.header.get("BSCALE", 1)
            + hdu.header.get("BZERO", 0)
        )
        if (maxerr := abs(check - arr).max()) > max_abs_err:
            raise ValueError(f"maxerr {maxerr}, above > {max_abs_err} cutoff")
    return hdu


def decimal_to_hourstring(min_ltst, max_ltst):
    return "_".join(gmap(lambda i: str(round(i * 24)), (min_ltst, max_ltst)))
