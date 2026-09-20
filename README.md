# layered-earth

実在DEMを5m標高帯に分け、連結成分・面積・接触グラフを解析する地域プロトタイプ。

**Status:** 解析・描画コード実装済み。実データ未取得、全球解析未対応。

## Run

```bash
python -m pip install -e '.[geo,test]'
python -m unittest discover -s tests -v
layered-earth --dem data/region.tif --source 'DATASET_URL' \
  --band-height 5 --connectivity 4 --output outputs/region
```

`outputs/region/results.md` にPNG/SVGの結果画像、CSVへのリンクを自動生成する。
実データ取得後、検証済み出力を `docs/results/` にコピーしてコミットできる。
元DEMはGit管理しない。出典・パラメータ・SHA-256は `provenance.json` に残す。

GeoTIFF: 単一バンド、非回転EPSG:4326、標高単位m、経度は昇順。
NPZ: `h` (行×列)、`lon_edges` (列+1)、`lat_edges` (行+1)、任意のboolean `mask`。
海抜0m以上という初期フィルタは陸海判定ではない。負標高も解析する場合は
`--min-elevation -500` 等を明示する。

[数理モデル・限界・研究計画](docs/model.md)
