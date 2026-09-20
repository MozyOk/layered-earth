"""Reproduce the committed Kanto experiment from the unmodified NOAA 60s TIFF.

Run from repository root: python scripts/kanto_experiment.py --dem data/ETOPO_2022_v1_60s_N90W180_surface.tif
Pillow fallback decodes the full raster (~1 GB); allow ~2 GB RAM.
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from PIL import Image, __version__ as pillow_version
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from layered_earth import analyze, write_csv, R_KM

SOURCE = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif'
SHA = '9d27d4b8ea8e76977e2988bca667d7c8fa68b927355feffcddd6b4875a7fd08e'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dem', type=Path, required=True)
    p.add_argument('--output', type=Path, default=Path('docs/results/kanto'))
    args = p.parse_args()
    digest = hashlib.sha256()
    with args.dem.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            digest.update(block)
    if digest.hexdigest() != SHA:
        raise ValueError('Unexpected input SHA-256; this experiment pins the NOAA raster')
    Image.MAX_IMAGE_PIXELS = None  # Only the checksum-pinned official file.
    with Image.open(args.dem) as im:
        if im.size != (21600,10800):
            raise ValueError('Unexpected raster dimensions')
        # Exact original cells, no interpolation: 138..141 E, 34..37 N.
        h = np.asarray(im.crop((19080,3180,19260,3360)), dtype=float)
    h[h == -99999] = np.nan
    lon,lat = np.linspace(138,141,181),np.linspace(37,34,181)
    valid = np.isfinite(h)&(h>=0)
    out=args.output
    out.mkdir(parents=True,exist_ok=True)
    comparisons=[]
    results={}
    for connectivity in (4,8):
        result=analyze(h,lon,lat,5,connectivity,valid)
        layers,labels,components,stats,edges=result
        results[connectivity]=result
        write_csv(out/f'layers-{connectivity}.csv',stats,['layer','min_m','components','area_km2'])
        write_csv(out/f'components-{connectivity}.csv',sorted(components,key=lambda c:-c['area_km2']),['id','layer','min_m','max_m','cells','area_km2'])
        write_csv(out/f'contacts-{connectivity}.csv',edges,['lower','upper','layer_jump'])
        # Independent per-row eligible cell count, not component aggregation.
        reference_area=sum(int(valid[r].sum())*R_KM**2*np.deg2rad(1/60)*
            abs(np.sin(np.deg2rad(lat[r]))-np.sin(np.deg2rad(lat[r+1]))) for r in range(180))
        assert np.isclose(sum(c['area_km2'] for c in components),reference_area)
        assert sum(c['cells'] for c in components)==int(valid.sum())
        assert all(e['layer_jump']>0 for e in edges)
        comparisons.append(dict(connectivity=connectivity,components=len(components),
            occupied_layers=len(stats),area_km2=reference_area,contacts=len(edges),
            singleton_components=sum(c['cells']==1 for c in components)))
    assert comparisons[1]['components']<=comparisons[0]['components']
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','svg.hashsalt':'layered-earth-kanto-v1'})
    def save(fig,name):
        fig.savefig(out/f'{name}.svg',metadata={'Date':None})
        fig.savefig(out/f'{name}.png',dpi=170)
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    for ax,k in zip(axes.flat,[0,1,2,4]):
        field=np.where(valid,1,0)
        field[valid & (results[4][0]==k)]=2
        ax.imshow(field,extent=[138,141,34,37],origin='upper',interpolation='nearest',
            cmap=ListedColormap(['#ffffff','#d9dde1','#2364aa']),vmin=0,vmax=2)
        ax.set_aspect(1/np.cos(np.deg2rad(35.5)))
        ax.set(title=f'{5*k} to <{5*(k+1)} m',xlabel='Longitude (E)',ylabel='Latitude (N)')
    fig.suptitle('Kanto | selected 5 m elevation bands\nETOPO 2022, 60 arc-sec; blue = selected, gray = other h >= 0, white = excluded',fontsize=12)
    save(fig,'bands')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for n,color,style in [(4,'#2364aa','-'),(8,'#b77520','--')]:
        stats=results[n][3]
        axes[0].plot([s['min_m'] for s in stats if s['min_m']<200],
                     [s['components'] for s in stats if s['min_m']<200],style,color=color,label=f'{n}-neighbor')
    axes[0].set(title='Component counts below 200 m',xlabel='Band lower edge (m)',ylabel='Components',ylim=(0,None))
    axes[0].legend()
    top=sorted(results[4][2],key=lambda c:-c['area_km2'])[:6][::-1]
    axes[1].barh([f"ID {c['id']} | {c['min_m']:g}–{c['max_m']:g} m" for c in top],
                 [c['area_km2'] for c in top],color='#2364aa')
    axes[1].set(title='Largest 6 components | 4-neighbor',xlabel='Spherical cell area (km²)')
    fig.suptitle('Kanto | DEM-dependent components, not physical islands')
    save(fig,'statistics')
    meta=dict(source=SOURCE,source_sha256=SHA,source_bytes=args.dem.stat().st_size,
        bounds_wsen=[138,34,141,37],pixel_window_xyxy=[19080,3180,19260,3360],shape=list(h.shape),
        horizontal_resolution_arcsec=60,vertical_datum='EGM2008',band_height_m=5,
        filter='finite elevation >= 0; no independent land mask',eligible_cells=int(valid.sum()),
        resampling='none',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        matplotlib=matplotlib.__version__,pillow=pillow_version),comparisons=comparisons)
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    rows='\n'.join(f"| {c['connectivity']} | {c['components']:,} | {c['occupied_layers']} | {c['area_km2']:,.2f} | {c['contacts']:,} | {c['singleton_components']:,} |" for c in comparisons)
    (out/'README.md').write_text(f'''# 関東周辺を5m標高帯に分ける — 初回実データ実験

ETOPO 2022 Surface / 60秒角。東経138–141°、北緯34–37°の180×180セルを無補間で抽出。
行政上の関東地方とは一致しない矩形領域。有効条件は標高0m以上、{valid.sum():,}セル。

![標高帯別の地図](bands.svg)

青が指定標高帯、灰色がそれ以外の標高0m以上、白が除外セル。
Plate Carrée座標で表示し、縦横比は北緯35.5°の距離比に合わせた。

## 計算結果

| 近傍 | 成分数 | 占有層数 | 対象面積 km² | 接触辺数 | 1セル成分数 |
|---|---:|---:|---:|---:|---:|
{rows}

![近傍比較と面積ランキング](statistics.svg)

対象面積は同じで、近傍規則を変えると成分数が変わる。飛地数はDEMと接続定義に依存する。
グラフは200m未満の層のみ表示。CSVにはすべての占有層を収録。

## 限界

- 水平解像度は60秒角（この緯度で概ね南北1.85km・東西1.51km）。5mの垂直精度ではない。
- 標高0m以上フィルタは陸海判定ではない。海面下の陸地を除外し、0m海セルを含み得る。
- 成分は領域境界で切断される。小さい成分・1セル成分を現実の独立した飛地と解釈しない。
- 接触辺は辺共有の地形接触。移動可能性・途中の未観測層を表さない。
- 全球解析・独立陸海マスク・垂直誤差モデルは未対応。

## 再現

```bash
python -m pip install -e . Pillow
mkdir -p data
curl -fL '{SOURCE}' -o data/ETOPO_2022_v1_60s_N90W180_surface.tif
python scripts/kanto_experiment.py --dem data/ETOPO_2022_v1_60s_N90W180_surface.tif
python -m unittest discover -s tests -v
```

約2GBのメモリを用意。入力SHA-256を検証してから処理し、不一致なら停止する。
PNG/SVG両方を生成するが、GitにはSVGを収録。元DEMはGit管理しない。

## 出典・検証

[NOAA配布元]({SOURCE}) / [ETOPO製品説明](https://www.ncei.noaa.gov/products/etopo-global-relief-model)

[パラメータ・SHA-256・実行環境](provenance.json)

- [4近傍の層別統計](layers-4.csv) / [面積ランキング](components-4.csv) / [接触グラフ](contacts-4.csv)
- [8近傍の層別統計](layers-8.csv) / [面積ランキング](components-8.csv) / [接触グラフ](contacts-8.csv)

成分のセル数合計、行別に独立計算した球面面積との一致、全接触辺の上向き性、
8近傍の成分数が4近傍以下になることを実データで検証した。
''',encoding='utf-8')
    print(json.dumps(meta,indent=2))


if __name__=='__main__':
    main()
