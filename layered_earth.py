"""Regional MVP. Explicit cell edges, spherical areas, no invented land mask."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

R_KM = 6371.0088


def analyze(h, lon_edges, lat_edges, band=5., connectivity=4, mask=None):
    h = np.asarray(h, dtype=float)
    lon, lat = np.asarray(lon_edges), np.asarray(lat_edges)
    if h.ndim != 2 or lon.shape != (h.shape[1]+1,) or lat.shape != (h.shape[0]+1,):
        raise ValueError("Expected h[rows,cols] and cell-edge coordinates")
    if not np.isfinite(band) or band <= 0 or connectivity not in (4, 8):
        raise ValueError("band must be positive; connectivity must be 4 or 8")
    if not (np.isfinite(lon).all() and np.isfinite(lat).all()):
        raise ValueError("Nonfinite coordinates")
    if not (np.all(np.diff(lon)>0) and (np.all(np.diff(lat)>0) or np.all(np.diff(lat)<0))):
        raise ValueError("Coordinates must be monotonic; longitude increasing")
    if np.max(abs(lat)) > 90 or lon[-1]-lon[0] >= 360-1e-8:
        raise ValueError("Regional MVP only: latitude within +/-90, longitude span <360")
    valid = np.isfinite(h)
    if mask is not None:
        if np.shape(mask) != h.shape:
            raise ValueError("Mask shape differs from DEM")
        valid &= np.asarray(mask, dtype=bool)
    if not valid.any():
        raise ValueError("No eligible cells")
    # NaN remains unassigned; negative elevation is NOT clipped to sea level.
    layers = np.full(h.shape, np.nan)
    layers[valid] = np.floor(h[valid]/band)
    area = R_KM**2 * np.abs(np.diff(np.sin(np.deg2rad(lat))))[:,None] * np.deg2rad(np.diff(lon))[None,:]
    labels = np.zeros(h.shape, dtype=np.int64)
    components, stats = [], []
    structure = ndimage.generate_binary_structure(2, 1 if connectivity == 4 else 2)
    offset = 0
    for k in np.unique(layers[valid]).astype(int):
        selected = valid & (layers == k)
        local, n = ndimage.label(selected, structure)
        labels[selected] = local[selected] + offset
        sizes = np.bincount(local.ravel(), weights=area.ravel(), minlength=n+1)
        counts = np.bincount(local.ravel(), minlength=n+1)
        for j in range(1, n+1):
            components.append(dict(id=offset+j, layer=int(k), min_m=k*band,
                                   max_m=(k+1)*band, cells=int(counts[j]), area_km2=float(sizes[j])))
        stats.append(dict(layer=int(k), min_m=k*band, components=n, area_km2=float(sizes[1:].sum())))
        offset += n
    # Edge-sharing contacts only. Corner contact does not define a shared boundary.
    pairs = []
    for a,b in [(labels[:,:-1],labels[:,1:]), (labels[:-1,:],labels[1:,:])]:
        keep = (a>0)&(b>0)&(a!=b)
        pairs.extend(zip(a[keep].tolist(),b[keep].tolist()))
    contacts = sorted({tuple(sorted(p)) for p in pairs})
    layer_by_id = {c['id']:c['layer'] for c in components}
    edges = []
    for a,b in contacts:
        if layer_by_id[a] > layer_by_id[b]:
            a,b = b,a
        edges.append(dict(lower=a, upper=b, layer_jump=layer_by_id[b]-layer_by_id[a]))
    return layers, labels, components, stats, edges


def read_dem(path):
    if path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as d:
            return d['h'], d['lon_edges'], d['lat_edges'], d['mask'] if 'mask' in d else None
    import rasterio
    with rasterio.open(path) as d:
        t = d.transform
        if d.crs != rasterio.crs.CRS.from_epsg(4326) or t.b or t.d or t.a <= 0:
            raise ValueError('GeoTIFF must be unrotated EPSG:4326, west to east')
        return (d.read(1, masked=True).astype(float).filled(np.nan),
                t.c+np.arange(d.width+1)*t.a, t.f+np.arange(d.height+1)*t.e, None)


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dem', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--source', required=True, help='Original dataset URL / provenance')
    p.add_argument('--band-height', type=float, default=5)
    p.add_argument('--connectivity', type=int, choices=[4,8], default=4)
    p.add_argument('--min-elevation', type=float, default=0,
                   help='Default excludes negative elevations; NOT a land mask')
    args = p.parse_args()
    h,lon,lat,mask = read_dem(args.dem)
    eligible = np.isfinite(h)&(h>=args.min_elevation)
    if mask is not None:
        eligible &= mask.astype(bool)
    layers,labels,components,stats,edges = analyze(h,lon,lat,args.band_height,args.connectivity,eligible)
    out = args.output
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'components.csv', sorted(components,key=lambda c:-c['area_km2']), ['id','layer','min_m','max_m','cells','area_km2'])
    write_csv(out/'layers.csv', stats, ['layer','min_m','components','area_km2'])
    write_csv(out/'contacts.csv', edges, ['lower','upper','layer_jump'])
    digest = hashlib.sha256()
    with args.dem.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            digest.update(block)
    meta = dict(source=args.source, sha256=digest.hexdigest(), shape=list(h.shape),
                band_height_m=args.band_height, connectivity=args.connectivity,
                min_elevation_m=args.min_elevation, land_mask_provided=mask is not None,
                area_model='sphere R=6371.0088 km', region_clipped=True,
                eligible_cells=int(eligible.sum()), components=len(components), contacts=len(edges))
    (out/'provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.fonttype':'none','font.size':11})
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    im=axes[0].pcolormesh(lon,lat,layers*args.band_height,cmap='cividis',rasterized=True)
    axes[0].set(xlabel='Longitude (degrees)',ylabel='Latitude (degrees)', title='Elevation-band lower edge')
    fig.colorbar(im,ax=axes[0],label='m')
    axes[1].plot([s['min_m'] for s in stats],[s['components'] for s in stats],color='#2364aa',marker='.',lw=1)
    axes[1].set(xlabel='Occupied band lower edge (m)',ylabel='Connected components',title=f'{args.connectivity}-neighbor connectivity',ylim=(0,None))
    axes[1].grid(alpha=.2)
    fig.suptitle('Layered Earth | regional DEM experiment')
    for extension in ['png','svg']:
        fig.savefig(out/f'overview.{extension}',dpi=180)
    plt.close(fig)
    (out/'results.md').write_text('# Regional results\n\n![DEM bands and component counts](overview.svg)\n\n'
        f'Components: {len(components):,}; eligible cells: {eligible.sum():,}.\n\n'
        'Source and parameters: [provenance](provenance.json). '
        'Region-clipped components; elevation filtering is not a land mask. '
        '5 m classes do not imply 5 m vertical accuracy.\n\n'
        '[Area ranking](components.csv) · [Layer statistics](layers.csv) · [Contact graph](contacts.csv)\n',encoding='utf-8')
    print(json.dumps(meta,indent=2))


if __name__ == '__main__':
    main()
