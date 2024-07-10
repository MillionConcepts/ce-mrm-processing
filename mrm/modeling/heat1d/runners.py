from itertools import product
import math
from multiprocessing import Pool
from pathlib import Path
import pickle
from sys import stdout
import time
from types import MappingProxyType as MPt
from typing import Mapping, Sequence

from dustgoggles.dynamic import exc_report
from hostess.monitors import Stopwatch

from mrm.modeling import heat1d


def model_metadata_dict(model: heat1d.Model):
    return {
        k: getattr(model, k)
        for k in (
            'albedo', 'lat', 'h', 'layer_scale_factor', 'nskin', 'nskinbot',
            'nyearseq', 'nperday', 'chi', 'fourier_mesh_number', 'dtsurf',
            'ztscale', 'eccentricity', 'obliquity', 'qb', 'ks', 'kd'
        )
    }


def model_data_dict(model: heat1d.Model):
    return {
        'temp': model.output_temp[:, 0:-1].T,
        'z': model.z[0:-1].T,
        'dz': model.dz.T,
        'rho': model.rho[0:-1].T,
        'ltst': model.lt
    }


def print_inline(*args, blanks=60):
    stdout.write(f"{' ' * blanks}\r")
    stdout.write(f"{' '.join(map(str, args))}\r")
    stdout.flush()


def run_model(lat, albedo, h, **params):
    model = heat1d.Model(lat=lat, albedo=albedo, h=h, **params)
    # print(f"----running {np.degrees(lat), albedo, h}----")
    model.run()
    # print(f"----completed {np.degrees(lat), albedo, h}----")
    return {'data': model_data_dict(model), 'meta': model_metadata_dict(model)}


def run_parallel(
    params: Mapping = MPt({}),
    latitudes: Sequence[float] = (0,),
    h_parameters: Sequence[float] = (0.08,),
    albedos: Sequence[float] = (0.06,),
    skip_existing: bool = False,
    threads: int = 4,
    outpath: Path = Path(".")
):
    start = time.time()
    results = {}
    pool = Pool(threads)
    points = tuple(product(latitudes, albedos, h_parameters))
    watch = Stopwatch()
    watch.start()
    for lat, alb, h in points:
        key = "_".join(map(str, (lat, alb, h))).replace('.', '')
        if (outpath / f"{key}.pkl").exists() and skip_existing is True:
            continue
        results[key] = (
            pool.apply_async(run_model, (math.radians(lat), alb, h), params)
        )
    n_runs = len(results)
    pool.close()
    Path("mat_output").mkdir(exist_ok=True)
    n_done, n_crashed = 0, 0
    while len(results) > 0:
        ready = [lah for lah, r in results.items() if r.ready()]
        for lah in ready:
            try:
                output = results[lah].get()
                with (outpath / f"{lah}.pkl").open("wb") as stream:
                    pickle.dump(output, stream)
                n_done += 1
            except KeyboardInterrupt:
                raise
            except Exception as ex:
                with (outpath / f"{lah}.err").open("wb") as stream:
                    pickle.dump(exc_report(ex), stream)
                n_crashed += 1
            finally:
                del results[lah]
        watch.update()
        n_complete, runtime = n_crashed + n_done, round(watch.total, 2)
        if n_complete > 0:
            s_per_item = runtime / n_done
            remaining = round(s_per_item * (n_runs - n_complete), 2)
        else:
            remaining = "NaN"
        print_inline(
            f"*** {n_crashed + n_done} / {n_runs} ({n_crashed} crashed); "
            f"{runtime}s (~{remaining}s remaining) ***"
        )
        time.sleep(1)
    pool.terminate()
    print(f"----done----")
    print(f"{time.time() - start}s")
    return results
