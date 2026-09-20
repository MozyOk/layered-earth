# 関東周辺を5m標高帯に分ける — 初回実データ実験

ETOPO 2022 Surface / 60秒角。東経138–141°、北緯34–37°の180×180セルを無補間で抽出。
行政上の関東地方とは一致しない矩形領域。有効条件は標高0m以上、18,350セル。

![標高帯別の地図](bands.svg)

青が指定標高帯、灰色がそれ以外の標高0m以上、白が除外セル。
Plate Carrée座標で表示し、縦横比は北緯35.5°の距離比に合わせた。

## 計算結果

| 近傍 | 成分数 | 占有層数 | 対象面積 km² | 接触辺数 | 1セル成分数 |
|---|---:|---:|---:|---:|---:|
| 4 | 14,563 | 511 | 50,970.97 | 29,663 | 13,345 |
| 8 | 13,475 | 511 | 50,970.97 | 27,542 | 11,941 |

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
curl -fL 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif' -o data/ETOPO_2022_v1_60s_N90W180_surface.tif
python scripts/kanto_experiment.py --dem data/ETOPO_2022_v1_60s_N90W180_surface.tif
python -m unittest discover -s tests -v
```

約2GBのメモリを用意。入力SHA-256を検証してから処理し、不一致なら停止する。
PNG/SVG両方を生成するが、GitにはSVGを収録。元DEMはGit管理しない。

## 出典・検証

[NOAA配布元](https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif) / [ETOPO製品説明](https://www.ncei.noaa.gov/products/etopo-global-relief-model)

[パラメータ・SHA-256・実行環境](provenance.json)

- [4近傍の層別統計](layers-4.csv) / [面積ランキング](components-4.csv) / [接触グラフ](contacts-4.csv)
- [8近傍の層別統計](layers-8.csv) / [面積ランキング](components-8.csv) / [接触グラフ](contacts-8.csv)

成分のセル数合計、行別に独立計算した球面面積との一致、全接触辺の上向き性、
8近傍の成分数が4近傍以下になることを実データで検証した。
