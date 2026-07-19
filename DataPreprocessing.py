# Antarctica/Greenland
import h5py
import numpy as np
from pyproj import CRS, Transformer
from pyresample import geometry as geom
from pyresample import kd_tree as kdt
from osgeo import gdal, osr
import os, glob
from scipy.spatial import cKDTree

def h5_tif(data_path, geo_path, out_folder):
    # Load reflectance data from the HDF file
    f = h5py.File(data_path, 'r')
    raw_ref = f['/Data/EV_1KM_RefSB'][:]
    cal_coef = f['/Calibration/VIS_Cal_Coeff'][:]
    esd = f.attrs['EarthSun Distance Ratio']

    # GEO files
    g = h5py.File(geo_path, 'r')
    lon = g['/Geolocation/Longitude'][:].astype(float)
    lat = g['/Geolocation/Latitude'][:].astype(float)
    solz = g['/Geolocation/SolarZenith'][:] * 0.01 * np.pi / 180.0
    # # Read additional angular data (in degrees)
    saa = g['/Geolocation/SolarAzimuth'][:] * 0.01
    vza = g['/Geolocation/SensorZenith'][:] * 0.01
    vaa = g['/Geolocation/SensorAzimuth'][:] * 0.01

    # # Calculate the Relative Azimuth Angle (RAA) = SAA − VAA
    raa = saa - vaa
    # Normalize RAA to the range [0, 180]
    raa = np.abs(raa)
    raa = np.where(raa > 180, 360.0 - raa, raa)

    # 
    lon = np.where((lon < -180) | (lon > 180), np.nan, lon)
    lat = np.where((lat < -90) | (lat > 90), np.nan, lat)


    # Antarctica
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3031", always_xy=True)
    # Greenland
    # transformer = Transformer.from_crs("EPSG:4326", "EPSG:3995", always_xy=True)
    x, y = transformer.transform(lon, lat)

    # Set the output extent and resolution
    llx, lly = np.nanmin(x), np.nanmin(y)
    urx, ury = np.nanmax(x), np.nanmax(y)
    ps = 1000
    cols = int(round((urx - llx) / ps))
    rows = int(round((ury - lly) / ps))
    areaExtent_m = (llx, lly, urx, ury)

    # Define gt  for unified use in subsequent file writes
    gt = [llx, ps, 0, ury, 0, -ps]


    # Regional definition
    areaDef = geom.AreaDefinition(
        'EPSG:3031', 'Antarctic Polar Stereographic', 'stere',
        CRS.from_epsg(3031), cols, rows, areaExtent_m
    )
    # areaDef = geom.AreaDefinition(
    #     'EPSG:3995', 'Arctic Polar Stereographic', 'stere',
    #     CRS.from_epsg(3995), cols, rows, areaExtent_m
    # )
    swath_def = geom.SwathDefinition(lons=lon, lats=lat)

    # Output grid points
    xs = llx + (np.arange(cols) + 0.5) * ps
    ys = ury - (np.arange(rows) + 0.5) * ps
    xv, yv = np.meshgrid(xs, ys)
    tgt_pts = np.column_stack((xv.ravel(), yv.ravel()))

    # KDTree
    valid_mask = np.isfinite(lon) & np.isfinite(lat)
    tree = cKDTree(np.column_stack((x[valid_mask], y[valid_mask])))

    # Radiometric calibration; six bands
    bands = {'B08': 3,  'B11': 6,'B13': 8,'B14':9, 'B15': 10, 'B19': 14}
    for bname, idx in bands.items():
        raw = raw_ref[idx, :, :].astype(float)
        if bname == 'B08':
            cal_idx =7        
        elif bname == 'B11':
            cal_idx = 10        
        elif bname == 'B13':
            cal_idx = 12
        elif bname == 'B14':
            cal_idx = 13
        elif bname == 'B15':
            cal_idx = 14
        elif bname == 'B19':
            cal_idx = 18        
        else:
            raise ValueError(f"Unknown band: {bname}")
    

        cos_solz = np.cos(solz)
        cos_solz[cos_solz == 0] = np.nan
        ref = (raw * cal_coef[cal_idx, 1] + cal_coef[cal_idx, 0]) / 100.0 * (esd ** 2) / cos_solz

    #  Print TOA max/min values 
        ref_min = np.nanmin(ref)
        ref_max = np.nanmax(ref)
        print(f"{bname} TOA reflectance: min={ref_min:.4f}, max={ref_max:.4f}")
        

        # Nearest-neighbor resampling
        result = kdt.resample_nearest(swath_def, ref, areaDef, radius_of_influence=5000, fill_value=np.nan)

        # fill NaN
        if np.isnan(result).sum() > 0:
            sval = ref[valid_mask]
            dist, idx_nn = tree.query(tgt_pts, k=1)
            filled_vals = sval[idx_nn].reshape((rows, cols))
            filled_vals[dist.reshape((rows, cols)) > 5000] = np.nan
            result = np.where(np.isnan(result), filled_vals, result)

        # Write  GeoTIFF
        gt = [llx, ps, 0, ury, 0, -ps]
        out_name = os.path.join(out_folder, os.path.basename(data_path).replace('.HDF', f'_{bname}_TOA_1km.tif'))
        d = gdal.GetDriverByName('GTiff').Create(out_name, cols, rows, 1, gdal.GDT_Float32)
        d.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3031)
        # srs.ImportFromEPSG(3995)
        d.SetProjection(srs.ExportToWkt())
        band_out = d.GetRasterBand(1)
        band_out.WriteArray(np.nan_to_num(result, nan=0.0))
        band_out.SetNoDataValue(-999)
        band_out.FlushCache()
        d = None

    # SZA / VZA
    for name, data in zip(['SZA', 'VZA'], [g['/Geolocation/SolarZenith'][:] * 0.01, g['/Geolocation/SensorZenith'][:] * 0.01]):
        data_res = kdt.resample_nearest(swath_def, data, areaDef, radius_of_influence=5000, fill_value=np.nan)
        if np.isnan(data_res).sum() > 0:
            sval = data[valid_mask]
            dist, idx_nn = tree.query(tgt_pts, k=1)
            filled_vals = sval[idx_nn].reshape((rows, cols))
            filled_vals[dist.reshape((rows, cols)) > 5000] = np.nan
            data_res = np.where(np.isnan(data_res), filled_vals, data_res)

        out_name = os.path.join(out_folder, os.path.basename(geo_path).replace('.HDF', f'_{name}_1km.tif'))
        d = gdal.GetDriverByName('GTiff').Create(out_name, cols, rows, 1, gdal.GDT_Float32)
        d.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3031)
        # srs.ImportFromEPSG(3995)
        d.SetProjection(srs.ExportToWkt())
        band_out = d.GetRasterBand(1)
        band_out.WriteArray(np.nan_to_num(data_res, nan=0.0))
        band_out.SetNoDataValue(-999)
        band_out.FlushCache()
        d = None
    print("TOA + SZA + VZA.tif generated.")
    #  SAA VAA RAA

    for name, data in zip(['SAA', 'VAA'], [g['/Geolocation/SolarAzimuth'][:] * 0.01, g['/Geolocation/SensorAzimuth'][:] * 0.01]):
        data_res = kdt.resample_nearest(swath_def, data, areaDef, radius_of_influence=5000, fill_value=np.nan)
        if np.isnan(data_res).sum() > 0:
            sval = data[valid_mask]
            dist, idx_nn = tree.query(tgt_pts, k=1)
            filled_vals = sval[idx_nn].reshape((rows, cols))
            filled_vals[dist.reshape((rows, cols)) > 5000] = np.nan
            data_res = np.where(np.isnan(data_res), filled_vals, data_res)

        out_name = os.path.join(out_folder, os.path.basename(geo_path).replace('.HDF', f'_{name}_1km.tif'))
        d = gdal.GetDriverByName('GTiff').Create(out_name, cols, rows, 1, gdal.GDT_Float32)
        d.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3031)
        # srs.ImportFromEPSG(3995)
        d.SetProjection(srs.ExportToWkt())
        band_out = d.GetRasterBand(1)
        band_out.WriteArray(np.nan_to_num(data_res, nan=0.0))
        band_out.SetNoDataValue(-999)
        band_out.FlushCache()
        d = None
    print("SAA + VAA.tif generated.")
    # ================= RAA (0–180) =================
    raa_res = kdt.resample_nearest(
        swath_def, raa, areaDef,
        radius_of_influence=5000, fill_value=np.nan
    )

    if np.isnan(raa_res).sum() > 0:
        sval = raa[valid_mask]
        dist, idx_nn = tree.query(tgt_pts, k=1)
        filled_vals = sval[idx_nn].reshape((rows, cols))
        filled_vals[dist.reshape((rows, cols)) > 5000] = np.nan
        raa_res = np.where(np.isnan(raa_res), filled_vals, raa_res)

    out_name = os.path.join(
        out_folder,
        os.path.basename(geo_path).replace('.HDF', '_RAA_1km.tif')
    )

    d = gdal.GetDriverByName('GTiff').Create(
        out_name, cols, rows, 1, gdal.GDT_Float32
    )
    d.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3031)
    # srs.ImportFromEPSG(3995)
    d.SetProjection(srs.ExportToWkt())

    band_out = d.GetRasterBand(1)
    band_out.WriteArray(np.nan_to_num(raa_res, nan=0.0))
    band_out.SetNoDataValue(-999)
    band_out.FlushCache()
    d = None

    print("RAA.tif generated.")        
    f.close()
    g.close()


if __name__ == "__main__":
    data_folder = r"\data"
    geo_folder = r"\geo"
    out_folder = geo_folder
    ms_files = sorted(glob.glob(os.path.join(data_folder, "*_1000M_MS.HDF")))
    for ms in ms_files:
        fname = os.path.basename(ms)
        geo = os.path.join(geo_folder, fname.replace("_1000M_MS.HDF", "_GEO1K_MS.HDF"))
        if os.path.exists(geo):
            print(f"\nprocessing: {fname}")
            h5_tif(ms, geo, out_folder)
        else:
            print(f"No corresponding GEO file found: {geo}")
