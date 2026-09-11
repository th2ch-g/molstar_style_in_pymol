# PyMOL で使う Mol* 風表示

[English](../guide.md) · [対応表と差分](../coverage.md) · [GPU/ray ギャラリー](../gallery.md)

`molstar_style` は独立した Python 実装です。Mol* の参照版は
`5b1b54ed03b03936041f514b33b8bb129b774d37`。Node.js、Mol*、CueMol、外部サーバーは不要です。
PyMOL 3.1 と OpenGL 2.1 / GLSL 1.20 に対応し、headless 環境でも ray 出力できます。

## 導入と基本操作

PyMOL が使う Python 環境へインストールします。

```sh
uv pip install --python <pymol-python> git+https://github.com/th2ch-g/molstar_style_in_pymol.git
```

PyMOL の Python または起動設定で登録します。mdtbx の `pymol_plugins` 経由なら自動登録されます。

```python
from molstar_style_in_pymol import __init_plugin__
__init_plugin__()
```

```text
molstar_style
molstar_style cartoon, selection=chain A, color=secondary-structure
molstar_style glossy, representation=ball-and-stick
molstar_style molecular-surface, transparency=0.4
molstar_style direct-volume, data=density.mrc
molstar_style segment, data=segments.npz
molstar_style interactions, selection=protein or ligand, name=contacts
molstar_style list
molstar_style refresh, name=all
molstar_style ray, filename=figure.png, width=1600, height=1200
molstar_style reset, name=all
```

既定の `polymer-and-ligand` は、ポリマーの cartoon、リガンド・イオンの球棒、水の半透明表示、
糖鎖の SNFG 記号を組み合わせます。`representation` で材質・効果に適用する形状を指定できます。
`color=keep` は元の原子色を使います。`state=0` は全状態を事前生成し、正の値は指定状態だけを表示します。
`quality` は lowest / lower / low / medium / high / higher / highest / auto / custom です。
`transparency=0..1` は各層の不透明度に乗算され、`keep` はプリセットの不透明度を保持します。

## 表示・色・設定

分子17種、ボリューム5種、粒子4種に加え、距離・角度・二面角、糖鎖、相互作用、架橋、膜、対称性、
DNATCO、衝突、軌道、トンネル、メッシュ、Kinemage、G3D、MVS、PAE パネルに対応します。
標準39色テーマと6サイズテーマ、拡張注釈用テーマの一覧は `molstar_style list` で確認できます。
形状別の `visuals` 名はパッケージの `reference.json` にあります。

`params` はローカル JSON ファイルまたは Python 辞書、`data` はローカル入力ファイルまたは辞書です。
PyMOL のコンマ区切り入力では JSON ファイルを使い、辞書は Python から渡してください。

```python
from molstar_style_in_pymol import molstar_style
molstar_style('cartoon', params={
    'tubularHelices': True,
    'material': {'metalness': 0.2, 'roughness': 0.4, 'bumpiness': 0},
    'postprocessing': {'occlusion': True, 'outline': True},
})
```

主な設定は、原子の `sizeFactor`・`multipleBonds`・`ignoreHydrogens`、cartoon の `aspectRatio`・
`arrowFactor`・`linearSegments`・`radialSegments`、表面の `resolution`・`probeRadius`・`smoothness`、
断面の `dimension`・`index`、密度の `isoValue`・`transferFunction`・`step` です。
色は `colorParams`、サイズは `sizeTheme` と `sizeParams` で設定します。
`clipPlanes` は最大6個の `[nx, ny, nz, offset]` を受け取り、内積と offset の和が非負の部分を残します。

## ローカル入力

CCP4/MRC、Cube、OpenDX、NPZ、CIF/BCIF、JSON、MVS JSON、G3D、Kinemage を読み込みます。
NPZ の密度は `values`、ボクセル番号から Å 座標への4×4行列は `transform` です。
読み込み済み PyMOL map は `selection=map_name` でも指定できます。
NPZ の pickle は無効です。ファイルや参照先の自動ダウンロードは行いません。

注釈の例:

```json
{"residues": [{"chain": "A", "resi": "10", "plddt": 94.5}]}
```

残基キーには `model`・`segi` を追加できます。重複する注釈はエラーにします。
原子配列は `atom_data` に置き、選択に対する `cmd.get_model` の原子順序・個数に合わせます。
mmCIF/BCIF の Model Archive 品質指標・entity 情報は自動対応付けします。
他の注釈形式は、上記の明示的な残基・原子スキーマへ変換して渡します。

| 表示 | 主な入力 |
| --- | --- |
| 距離・角度・二面角 | `positions` に2/3/4点、または選択原子の0始まり `indices` |
| 任意ラベル | `text`、`positions`。位置省略時は選択の中心 |
| メッシュ | `vertices`、三角形 `faces`、任意の `colors`・`opacity`・`transform` |
| 粒子 | `particles` の各行に `position`・`radius`。姿勢は XYZW の `quaternion` または `axes` |
| 繊維・標的 | 粒子行に `points` / `target` を追加 |
| 架橋・接触・衝突 | `cross_links` / `interactions` / `clashes` に原子 `indices`、種別、距離上下限など |
| 膜 | `membrane` に `center`・`normal`・`thickness`・`radius` |
| 対称性 | `symmetry.axes` に `start`・`end`・`order`、任意の cage 頂点と辺 |
| DNATCO | `steps` に `class`・`score`・`positions`。confal の点順は O3′、P、OP1、OP2、O5′ |
| トンネル | `positions` と正の `radii` |
| 軌道 | Mol* 形式の `basis.atoms`・`orbitals`、または計算済み Cube |
| PAE | 正方行列 `predicted_aligned_error` / `matrix`。欠損は null |
| 単位格子 | `cell: [a,b,c,alpha,beta,gamma]`、または PyMOL の対称性情報 |

相互作用は入力を省くと幾何的に推定します。明示的な水素があれば水素結合の角度条件も確認します。
ASA は Shrake–Rupley 法で計算します。衝突の推定は `compute=true` を指定した場合のみです。
膜・対称性・DNATCO・品質評価の科学的注釈は与えられた結果を可視化し、欠損値を推測しません。

軌道は実球面調和関数 L=0..4、gaussian / cca / cca-reverse 順に対応し、基底の中心は Bohr 単位です。
電子密度は占有数を掛けた軌道の二乗和です。G3D はローカル解像度ブロックを読み、染色体・ハプロタイプ・
領域で絞れます。Kinemage はベクトル、リボン、三角形、球、点、ラベルを読み込みます。

MVS はローカル構造の読み込み、成分選択、表現、色、不透明度、変換、ラベル、snapshot 選択に対応します。
複雑な注釈ノードは明示スキーマへ変換してください。未対応ノードを無言で省略することはありません。
Mol* のブラウザーアプリ、サーバー、編集 UI、全 MVS プロトコルを移植したものではありません。

PAE は Qt パネルに表示し、セルをクリックすると対応する2残基を選択します。
行列と選択残基数の一致が必要です。パネルのみのビューでは PNG/ray 操作で行列画像を書き出します。

## 描画の差分と復元

材質は matte / plastic / glossy / metallic。効果は occlusion / outline / shadow / cel / xray /
unlit / bloom / dof / illumination / antialias です。illumination は画面空間の間接光近似で、
Mol* の progressive path tracing そのものではありません。
背景は横グラデーション、放射状グラデーション、ローカル画像、6面の skybox を指定できます。
明示的に指定しない限り、既存の背景設定は変更しません。

不透明メッシュは GLSL、透明メッシュは頂点照明を焼き付けた CGO、密度は GPU ray marching で表示します。
標準 `ray` にも保持したメッシュと3方向の密度積分投影が含まれます。既定の PyMOL 透明度モードは重なる透明面を積算しないため、標準 ray はこの粗い投影を使い、専用 ray は積算用設定を一時適用して復元します。専用 `molstar_style ray` は視点に合わせた
密度断面、輪郭、背景合成、画面効果を追加します。native ray の効果はメッシュ深度サンプルによる近似です。
表面・blob・相互作用などの数値差もあるため、厳密な Mol* 画像の一致は保証しません。差分は対応表に記載しています。

同じ原子を複数の名前で表示できます。最後のビューを reset すると元の representation が戻ります。
同名の再適用では新しい形状の生成・読み込み完了後に置き換えます。失敗時は以前のビューを保持します。
元の座標・結合・色・二次構造は変更しません。編集後は `refresh` してください。
セッションには復元情報を保存し、再読込時に形状を再生成します。必要なローカルファイルが失われても元表示を復元します。

全状態を事前生成します。既定 `cache_mb=2048` は形状と native CGO の float データ量、
`gpu_cache_mb=256` は VBO・密度テクスチャを制限します。Python・PyMOL の管理領域、原子スナップショット、
画面サイズに応じた FBO は別途必要です。測定スクリプトは RSS も記録します。

```sh
pixi run test
pixi run check
pixi run uv run --no-project python tests/check_real.py --gui
pixi run uv run --no-project python tests/check_real.py --gui --visuals
pixi run uv run --no-project python tests/check_real.py --gui --benchmark --structure local_protein.cif
```

検証は別の PyMOL プロセスで行います。再生成可能な画像・性能結果は `.cache/` に出力され、Git 対象外です。

ASCII ラベルは同梱ベクトルフォントを使います。Unicode ラベルは Qt のフォントを使い、
headless の場合は必要な字形を持つローカル TTF/OTF を `params.font` で指定します。
文字形状も事前生成し、GPU と ray で共有します。
