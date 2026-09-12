# Mol*風のPyMOLスタイル

[Package README](../../README.md) | [English](../pymol_molstar.md)

`molstar_style_in_pymol` は分子形状、ボリューム、粒子、注釈、GPU 材質を追加する
独立した Python パッケージです。Python 3.10+、NumPy、SciPy、PyOpenGL、Pillow、
Gemmi、msgpack、scikit-image が必要です。対話表示には PyMOL 3.1 と Qt、
互換 OpenGL 2.1 / GLSL 1.20 コンテキストを使います。`pixi install` で PyMOL を含む
開発環境を導入できます。Node.js、Mol*、CueMol、mdtbx、外部サーバーは実行時に不要です。
PyMOL のソースコードと標準コマンドは変更しません。

## 使い始める

このリポジトリで `pixi install` を実行し、`pixi run pymol` で起動します。
既存の PyMOL 環境へ入れる場合は、その Python インタープリターを
`<pymol-python>` に指定してください。

```sh
uv pip install --python <pymol-python> "git+https://github.com/th2ch-g/molstar_style_in_pymol.git"
```

PyMOL の Python コンソール、または `.pymolrc.py` でコマンドを登録します。

```python
from molstar_style_in_pymol import __init_plugin__

__init_plugin__()
```

パッケージを import するだけでは PyMOL の import・起動・設定変更は行いません。
mdtbx の `pymol_plugins` 連携では自動登録されます。
構造を読み込んだ後、PyMOL コマンドラインで実行します。

```text
molstar_style
molstar_style cartoon, selection=chain A, color=secondary-structure
molstar_style glossy, representation=ball-and-stick
molstar_style molecular-surface, transparency=0.4
molstar_style list
molstar_style help
```

既定のビュー名は `molstar` です。同じ名前でスタイルを適用すると、そのビューを置き換えます。
既定の `polymer-and-ligand` は、ポリマーの cartoon、リガンド・イオンの球棒、
水の半透明表示、糖鎖の SNFG 記号を組み合わせます。`protein-and-nucleic` はポリマーだけを表示します。

## スタイルと設定

スタイル名

| 分類 | 名前 |
| --- | --- |
| 分子形状 | `cartoon`, `backbone`, `ball-and-stick`, `blob-surface`, `carbohydrate`, `ellipsoid`, `gaussian-surface`, `gaussian-volume`, `label`, `line`, `molecular-surface`, `orientation`, `plane`, `point`, `putty`, `spacefill`, `polyhedron` |
| ボリューム | `direct-volume`, `dot`, `isosurface`, `segment`, `slice` |
| 粒子 | `particle-spacefill`, `particle-orientation`, `particle-fibers`, `particle-target` |
| 測定と形状 | `distance`, `angle`, `dihedral`, `shape-label`, `shape-orientation`, `shape-plane`, `unitcell` |
| 拡張表示 | `interactions`, `cross-link-restraint`, `membrane-orientation`, `assembly-symmetry`, `confal-pyramids`, `ntc-tube`, `clashes`, `orbital`, `orbital-density`, `tunnel`, `mesh`, `kinemage`, `g3d`, `mvs`, `pairwise-metric`, `annotation-label`, `custom-label` |
| 材質 | `matte`, `plastic`, `glossy`, `metallic` |
| 画面効果 | `outline`, `occlusion`, `shadow`, `cel`, `xray`, `unlit`, `bloom`, `dof`, `illumination`, `background`, `antialias` |
| プリセット | `default`, `auto`, `empty`, `polymer-and-ligand`, `protein-and-nucleic`, `polymer-cartoon`, `atomic-detail`, `coarse-surface`, `illustrative`, `auto-lod`, `mesoscale`, `validation-geometry`, `validation-density`, `validation-rci`, `quality-plddt`, `quality-qmean`, `partial-charges` |

`nucleic` は `cartoon` の別名です。`empty` は管理対象の形状を消去します。
[ギャラリー](../gallery.md)に全表示と必要な入力を掲載しています。

```text
molstar_style [style], selection=all, representation=auto, color=auto, quality=medium, name=molstar
```

`representation` で材質・効果に適用する形状を指定します。`molstar_style list` は
スタイル、標準39色テーマ、6サイズテーマ、拡張注釈用テーマの一覧を表示します。
`quality` は lowest / lower / low / medium / high / higher / highest / auto / custom です。
独自の分割数や解像度は `params` で指定します。
`transparency=0..1` は各層の不透明度に乗算され、`keep` はプリセットの不透明度を保持します。
形状別の `visuals` 名はパッケージの `reference.json` にあります。
必要な科学的注釈が欠けている場合は、既存ビューを置き換える前にエラーにします。

### 表示パラメーター

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
`clipPlanes` は最大6個の `[nx, ny, nz, offset]` を受け取り、内積と offset の和が非負の部分を残します。

## 既定の配色

`color=auto` は分子の表示形式に応じて色テーマを選びます。原子表示では元素色を使い、
炭素は鎖色を引き継ぎます。ポリマーの cartoon と分子表面は鎖色を使います。

| モード | 配色 |
| --- | --- |
| `auto` | 表示形式に応じた分子色 |
| `keep` | 既存の PyMOL 原子色 |
| `chain-id` | 元オブジェクト全体の順序に基づく鎖色 |
| `secondary-structure` | ヘリックス・シート・コイルの色 |
| `element-symbol` | 元素色。炭素は既定で鎖色 |
| `uniform` | `colorParams.value` で指定する単色 |

```text
molstar_style cartoon, color=keep
molstar_style cartoon, color=chain-id
molstar_style cartoon, color=secondary-structure
```

色は `colorParams`、サイズは `sizeTheme` と `sizeParams` で設定します。
元の原子色・座標・物性は変更しません。原子色の編集後は `refresh` が必要で、
その原子色を使うテーマにのみ反映されます。

`chain-id` は部分選択でも元オブジェクト全体の配色順を保持します。mmCIF では PyMOL が保持する
entity ID (`custom`) と入力順 (`rank`) から、Mol* と同じ分子種ごとの鎖順を復元します。
PDB 入力から欠けている mmCIF の情報を復元することはできません。

注釈の例:

```json
{"residues": [{"chain": "A", "resi": "10", "plddt": 94.5}]}
```

残基キーには `model`・`segi` を追加できます。重複する注釈はエラーにします。
原子配列は `atom_data` に置き、選択に対する `cmd.get_model` の原子順序・個数に合わせます。
mmCIF/BCIF の Model Archive 品質指標・entity 情報は自動対応付けします。
他の注釈形式は、上記の明示的な残基・原子スキーマへ変換して渡します。

## ローカル入力

CCP4/MRC、Cube、OpenDX、NPZ、CIF/BCIF、JSON、MVS JSON、G3D、Kinemage を読み込みます。
NPZ の密度は `values`、ボクセル番号から Å 座標への4×4行列は `transform` です。
読み込み済み PyMOL map は `selection=map_name` でも指定できます。
NPZ の pickle は無効です。ファイルや参照先の自動ダウンロードは行いません。

| 表示 | 主な入力 |
| --- | --- |
| 距離・角度・二面角 | `positions` に2/3/4点、または選択原子の0始まり `indices` |
| 任意ラベル | `text`、`positions`。位置省略時は選択の中心 |
| メッシュ | `vertices`、三角形 `faces`、任意の `colors`・`opacity`・`transform` |
| 粒子 | `particles` の各行に `position`・`radius`。姿勢は XYZW の `quaternion` または `axes` |
| 繊維 | 粒子行に `points`。`linearSegments` と `tubeSizeFactor` で補間と太さを指定 |
| 粒子 target | `target` に `targets` 辞書のキーを指定。対象は shape / structure / volume |
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

軌道の既定 isovalue は最大絶対値の15%です。Mol* の累積確率による閾値計算とは異なります。
粒子の `scale` は3軸の正の無次元倍率で、spacefill では `radius * sizeFactor` に乗算します。
粒子 orientation の軸長は既定で10 Åです。target は指定形状を粒子ごとに回転・配置します。
旧実装の XYZ 終点への矢印は誤りのため廃止しました。

```json
{"targets": {"unit": {"kind": "shape", "vertices": [[-1,-1,0],[1,-1,0],[0,1,0]], "faces": [[0,1,2]]}},
 "particles": [{"target": "unit", "position": [5,0,0], "radius": 2, "quaternion": [0,0,0,1]}]}
```

shape はメッシュ入力、structure は原子 `positions`・`radii` と `type=spacefill/blob-surface`、
volume はローカル `grid` と `type=isosurface/dot` を使います。回転中心は `center`、省略時は境界箱の中心です。
`scaleByRadius` は shape で既定 true、他は false。`targetColor=source` で対象の色を保持します。
Mol* の target 配信・動的 LOD は含みません。[全表示の監査結果](../audit.md)に検証内容と残る近似を記載しています。

MVS はローカル構造の読み込み、成分選択、表現、色、不透明度、変換、ラベル、snapshot 選択に対応します。
複雑な注釈ノードは明示スキーマへ変換してください。未対応ノードを無言で省略することはありません。
Mol* のブラウザーアプリ、サーバー、編集 UI、全 MVS プロトコルを移植したものではありません。

PAE は Qt パネルに表示し、セルをクリックすると対応する2残基を選択します。
行列と選択残基数の一致が必要です。パネルのみのビューでは PNG/ray 操作で行列画像を書き出します。

```text
molstar_style isosurface, selection=density
molstar_style direct-volume, data=density.mrc
molstar_style segment, data=segments.npz
molstar_style interactions, selection=protein or ligand, name=contacts
```

上の例では map `density` または入力ファイルを事前に用意します。
相互作用の例は、既存の `protein`・`ligand` 選択を使います。

ASCII ラベルは同梱ベクトルフォントを使います。Unicode ラベルは Qt のフォントを使い、
headless の場合は必要な字形を持つローカル TTF/OTF を `params.font` で指定します。
文字形状も事前生成し、GPU と ray で共有します。

## 選択とstate管理

`selection` は PyMOL の分子選択または読み込み済み map を受け取ります。
`state=0` は全状態を事前生成し、正の値は指定状態だけを表示します。
複数のビューを保持するには別々の `name` を指定します。

```text
molstar_style refresh
molstar_style refresh, name=all
molstar_style reset
molstar_style reset, name=all
```

`name` を省略すると既定名 `molstar` のビューだけが対象になります。`name=all` は全ビューに適用します。
特定のビューには `molstar_style refresh, name=contacts` または
`molstar_style reset, name=contacts` のように名前を指定してください。

同じ原子を複数の名前で表示できます。最後のビューを reset すると元の representation が戻ります。
同名の再適用では新しい形状の生成・読み込み完了後に置き換えます。失敗時は以前のビューを保持します。
元の座標・結合・色・二次構造は変更しません。編集後は `refresh` してください。
セッションには復元情報を保存し、再読込時に形状を再生成します。必要なローカルファイルが失われても元表示を復元します。

`state=0` では全状態を事前生成します。既定 `cache_mb=2048` は形状と native CGO の float データ量、
`gpu_cache_mb=256` は VBO・密度テクスチャを制限します。Python・PyMOL の管理領域、原子スナップショット、
画面サイズに応じた FBO は別途必要です。測定スクリプトは RSS も記録します。

## 画像出力

```text
ray 1600, 1200
png figure_ray.png
molstar_style png, filename=figure.png, width=1600, height=1200
molstar_style ray, filename=figure_ray.png, width=1600, height=1200
```

専用 `png` 操作は Qt GUI の GPU 表示を出力します。
標準 ray と専用 ray は headless PyMOL でも使えます。

標準 `ray` にも保持したメッシュと3方向の密度積分投影が含まれます。
既定の PyMOL 透明度モードは重なる透明面を積算しないため、標準 ray はこの粗い投影を使い、
専用 ray は積算用設定を一時適用して復元します。専用 `molstar_style ray` は視点に合わせた
密度断面、輪郭、背景合成、画面効果を追加します。
管理対象の不透明メッシュだけを含むシーンでは、専用 ray は各画素の密度を直接積算し、メッシュ深度で遮蔽します。
管理対象外の形状や透明メッシュを含む場合は native ray の断面近似を使います。
表面・blob・相互作用などの数値差もあるため、厳密な Mol* 画像の一致は保証しません。差分は対応表に記載しています。
専用 ray は現在の視点で同じ材質計算を各頂点に適用します。管理対象の形状だけを表示するシーンでは
PyMOL の追加照明を一時的に無効化して二重照明を防ぎ、設定を全て復元します。
管理対象外の形状が見える混在シーンでは、その照明を維持するため管理対象の色も影響を受けます。
標準 `ray` は保持済み頂点色と PyMOL の照明を使うため、視点変更後は `refresh` が必要です。
専用 ray はフラグメント微分による粗さ補正と bump を省略します。輪郭の平滑化・頂点照明・透明度の
サンプリングにも差が残ります。[比較画像と数値](fidelity.md)を参照してください。

## 表示形式の照合結果

参照は Mol* の固定コミット `5b1b54ed03b03936041f514b33b8bb129b774d37` です。
[形状の監査](../audit.md)に修正内容、[対応表](../coverage.md)に各表示の入力と制約を記載しています。
[実際の Mol* との比較](fidelity.md)では入力とカメラを揃えて検証しています。

| Mol* 参照画像（8GNG） | PyMOL GPU | PyMOL 専用 ray |
| --- | --- | --- |
| ![Molstar 8GNG](../gallery/reference-molstar-8gng-color-gpu.png) | ![PyMOL GPU 8GNG](../gallery/reference-pymol-8gng-color-gpu.png) | ![PyMOL ray 8GNG](../gallery/reference-pymol-8gng-color-ray.png) |

楕円体は異方性変位テンソル（Å²）が必須です。未定義の原子は省略し、全原子で未定義ならエラーにします。
既定の半軸長は参照 Mol* と同じ `1.5958 * sqrt(abs(固有値)) * sizeFactor` です。
この定数は厳密な50%確率の係数とは異なります。`probability` を明示した場合のみ3次元χ²分布の分位点を使います。
等しい固有値ではテンソル由来の球になります。結合半径は既定で0.1 Åです。

backbone の既定半径は0.3 Å、cartoon はアスペクト比適用前で0.2 Åです。
cartoon は Mol* の残基ごとの Catmull–Rom 曲線、ペプチド方向からのフレーム、シート平滑化、
平らな β シートと矢尻を再現します。`helixProfile` / `nucleicProfile` は
`elliptical`・`rounded`・`square` を指定でき、既定はそれぞれ楕円・矩形です。
核酸の主鎖は O3' を使います。`tubularHelices` は helixorient 中心軸、`roundCap` はその丸い末端です。
[実際の Mol* との比較](fidelity.md)で曲線・表面頂点・画像の差を確認できます。
putty は既定で `0.2 * (0.2 + 0.1 * B)` Åを使い、`sizeTheme` も反映します。
旧平方根式は `bfactorScale` を明示した場合のみ使います。
構造の orientation は楕円体、shape-orientation は配向ボックスが既定です。
それぞれ `sizeFactor` / `scaleFactor` で大きさを変えられます。

構造の `plane` は原子を色分けした断面です。`imageResolution`、`axis=a/b/c`、`offset`、`margin`、
`cutout`、または `mode=plane` と `plane={point: [...], normal: [...]}` で指定します。
単純な当てはめ平面は `shape-plane` です。ボリュームの `slice` は格子の x/y/z 軸に対応し、既定は x 軸の中央です。
uniform 色では密度で明度を変え、`isoValue` より小さい値を透過させます。斜め断面・周期写像は未対応です。

線・点はワールド座標での近似です。線半径は `0.02 * sizeFactor * size`、点半径は `0.075 * sizeFactor * size` Åです。
線の既定 sizeFactor は2、点は1です。attenuation の真偽値を半径には使いません。
画面ピクセル寸法・カメラ距離による減衰・全ての上流結合フィルターは未実装です。
測定の `arcScale` は短い方の腕に対する割合で、角度・二面角は既定で扇形を描きます。
volume dot の既定半径は1 Åです。負の等値面では閾値以下の点を選びます。
Gaussian volume の分子色は GPU と ray の両方に反映します。

## 材質と輪郭線の照合

既定の材質係数は参照 Mol* に合わせています。`material` 辞書では
`metalness`・`roughness`・`bumpiness` を `[0,1]` で指定できます。

| 材質 | Metalness | Roughness | Bumpiness |
| --- | ---: | ---: | ---: |
| `matte` | 0.0 | 1.0 | 0.0 |
| `plastic` | 0.0 | 0.2 | 0.0 |
| `glossy` | 0.0 | 0.6 | 0.0 |
| `metallic` | 1.0 | 0.6 | 0.0 |

効果は occlusion / outline / shadow / cel / xray /
unlit / bloom / dof / illumination / antialias です。illumination は画面空間の間接光近似で、
Mol* の progressive path tracing そのものではありません。
背景は横グラデーション、放射状グラデーション、ローカル画像、6面の skybox を指定できます。
明示的に指定しない限り、既存の背景設定は変更しません。

不透明・透明メッシュとも Mol* の GGX / Schlick / Smith 材質計算を GLSL で行います。
透明な三角形は奥から描画し、密度は GPU ray marching で表示します。
既定の光源はカメラに対して inclination=150°、azimuth=320°、強度0.6、環境光は0.4です。
`params.lighting` で `light`（最大8灯、それぞれ角度・`color`・`intensity`）、
`ambientColor`・`ambientIntensity`・`exposure` を指定できます。
GPU は `flatShaded` と Mol* のノイズを使った `bumpFrequency` / `bumpAmplitude` にも対応します。
`bumpFrequency` の既定は0です。

occlusion は Mol* と同じ32個の固定サンプル、視点座標の深度法線、両側ぼかしを使います。
GPU と専用 ray で `radius`（Å の log2、既定5）、`bias`（0.8）、`blurKernelSize`（15）、
`blurDepthBias`（0.5 Å）を指定できます。多段解像度と透明物体の SSAO は未対応です。
積算密度には単一の表面深度がないため、直接ボリュームを含む層では表面 SSAO を省きます。

## 検証

```sh
pixi install --locked
pixi run test
pixi run check
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_real.py --gui --output .cache/validation
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_real.py --gui --visuals --output .cache/validation-visuals
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_real.py --gui --benchmark --structure local_protein.cif --output .cache/benchmark
uv build --python .pixi/envs/default/bin/python
```

検証スクリプトは GPU/ray 画像、ローカル入力、元表示の復元を確認します。
`--visuals` は標準の subvisual を個別に確認します。`--benchmark` はローカルの
タンパク質 CIF/PDB から500残基・100状態の合成データを作るため、`--structure` で入力を指定します。
GUI 検証は別の PyMOL プロセスで実行し、利用中のセッションは使いません。

Mol* と PyMOL の形状・材質・カメラを揃える検証は[比較結果](fidelity.md)と
[再実行手順](../../tests/reference/README.md)を参照してください。参照側のツールは検証専用です。
ギャラリー画像の再生成手順は[ギャラリー](../gallery.md#regenerate)にあります。
一時画像・環境依存の性能結果は `.cache/` に置き、Git 対象外とします。
公開するギャラリーと比較画像だけを `docs/gallery/` に保存し、CI を使わずローカルで描画します。
第三者のコード・データの帰属は [NOTICE](../../NOTICE) に記載しています。
