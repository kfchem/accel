# ACCeL v2.0 基盤監査レポート（v0.3.1 時点）

本書は ACCeL v2.0（`Structure` / `Flow` / 増分再計算）に向けて、現行コードベース
（v0.3.1, commit `d746482`）を監査した結果と、安全な移行のための設計論点をまとめたものです。
**本書は仕様の確定ではありません。** 科学的意味・既存ユーザー挙動・Flow UX・同一性判定に関わる
事項は「論点」として残し、明示的な決定を待ちます（§8）。

ファイル参照は `path:line` 形式（v0.3.1 の行番号）。

---

## 1. 現行アーキテクチャ

```text
accel/__init__.py            Box, System, Systems, Atom, Atoms を公開
accel/base/
  atoms.py      Atom / Atoms / Bonds / BondType     ← 原子・座標・結合・幾何計測
  systems.py    System / Systems                    ← 1構造 + 状態 + メタデータ / そのリスト
  boxcore.py    BoxCore                             ← chain API 本体（~40メソッド）
  box.py        Box(BoxCore)                        ← filetype 振り分けメソッド + .plugin
  selector.py   FuncSelector / Selectors            ← filetype → parser/実行関数 のレジストリ
  formats.py    xyz / mol 入出力, テンプレートエンジン
  topology.py   rmsdpruning, map_numbers
  modeler.py    Modeler（結合・芳香族・立体・対称行列・マッピング）, Chains（置換基順位付け）
  tools.py      change_dir, float_to_str
accel/plugin/
  gaulib.py (Gaussian)  orclib.py (ORCA)  xtblib.py (xTB)  maelib.py (Maestro)
  pbslib.py (PBS/qsub)  txtlib.py (テキスト編集)
accel/util/
  constants.py (Elements[pickle], CONSTANTS, Units)  datadict.py (Data)  filetype.py (FileType)
  log.py (グローバル logger, atexit で .mcl 出力)  execmd.py (実行コマンド名マップ)  matrix.py  dialog*.py
```

中心的な仕組み:

- **import 時登録**: 各 plugin モジュールの `@FileType.add(...)` と `@Selectors.<op>.add(filetype)` が
  import 時にグローバルレジストリへ登録される（`box.py` → `accel.plugin` → 全 plugin import）。
- **filetype 判定**: `System.path` 代入時に `FileType.analyse` が実ファイルを開いて判定（metric 昇順で最初に
  True を返した型）。`systems.py:62-105` の `__setattr__` 内で I/O が発生する。
- **振り分け**: `Box.read_energy()` 等は `adaptive_function_caller`（`box.py:11`）で active な System を
  filetype ごとに分け、`selector(ft, Box().bind(confs))` として plugin の**非束縛メソッド**を
  「素の Box」を self にして呼ぶ。該当関数がなければ `deactivate("could not find an appropriate function")`。
- **状態**: 全ての処理は `System` オブジェクトを in-place で書き換え、`self` を返す。

## 2. System / Systems の責務

`System`（`systems.py:14`）は `__slots__` で以下を持つ:

| 属性 | 性質 | 備考 |
|---|---|---|
| `atoms: Atoms` | 科学データ | 元素・座標・原子形式電荷・結合（`atoms.bonds`）・stereo |
| `total_charge` / `charge` | 科学データ | `charge` は `total_charge` が None なら原子形式電荷の和（`systems.py:107`） |
| `multiplicity` | 科学データ | **既定値 1**（読み込み元に情報がなくても 1 が入る） |
| `energy` | 計算結果 | kcal/mol の単一スカラー。SCF / 自由エネルギー / 相対値のいずれかは**記録されない** |
| `data: Data` | 計算結果＋任意メタデータ | `g16_scf`, 熱補正, 振動, `rotamer` 行列, `sdf` 生データ, `jobid` 等が混在 |
| `path`, `name`, `filetype` | 由来/ファイル | `name` は既定でファイル stem、`path` 代入で filetype 再判定 |
| `label` | グルーピング | 分子（配座集合）単位の処理キー |
| `state`, `history` | ワークフロー状態 | フィルタ結果（active/inactive）と理由の文字列履歴 |
| `distribution` | 集団依存の解析結果 | Boltzmann 分布（同じラベル集団に依存） |
| `cache` | 一時領域 | `relative_energy`, `distr_factor` 等 |

加えて `System.__setattr__` は型変換・DEBUG ログ・filetype 再判定・`Log.input_dir` 更新（グローバル）
を行う。同一性（identity）の概念はなく、Python オブジェクト同一性と `name` のみ。

`Systems`（`systems.py:185`）は `MutableSequence` で、`has_state/has_label/has_energy/has_data/
has_filetype/has_bonds/has_distribution` の**ビュー**（同じ System オブジェクトを共有する新 Systems）、
`labels`/`filetypes`（ソート済み dict）、`get()`（None=最小エネルギー or 名前順先頭、int=1始まり、
str=部分一致）、`sorted()`、`duplicate()` を提供する。

## 3. Box の責務

`BoxCore`（`boxcore.py:37`）は `contents: Systems` と `data` を持つ collection wrapper。責務の混在:

| 分類 | メソッド |
|---|---|
| 取り込み | `__init__`(glob/dir/path/System/Systems/Box), `add`, `bind` |
| 状態・ラベル・メタ | `get`, `labeling`, `set_label`, `set_state`, `set_data`, `zero_fill`, `count`, `show`, `duplicate` |
| エネルギー解析 | `energy_limit`, `calc_rel_energy`, `calc_distribution`, `calc_energy`, `only_minimum`, `get_average` |
| 構造アルゴリズム | `calc_bonds`, `calc_stereo`, `calc_symm`, `rmsd_limit`, `map_numbers`, `calc_length/angle/dihedral`, `modify_length`, `convert_to_mirror` |
| ファイル I/O | `read_xyz`, `write_xyz`, `read_mol`, `write_mol`, `write_input`, `export_data`, `copy_files`, `search` |
| filetype 振り分け（Box） | `read_atoms`, `read_energy`, `read_thermal`, `check_end`, `check_freq`, `calc_free_energy`, `run`, `submit` |
| plugin | `.plugin.gau/.orc/.xtb/.mae/.pbs/.txt`（同じ contents を共有する plugin Box を返す） |

chain を支える実装上の約束:

- **全メソッドが `self` を返す**（`get`/`get_average`/`get_resub`/`get_irc`/`get_trj`/`get_unzip` を除く）。
- **処理対象は原則 `self.get()`（active のみ）**。ただし `labeling`/`set_label`/`set_data`/`set_state` は
  inactive を含む全件に作用する（`boxcore.py:83-106`）。
- **フィルタは削除ではなく deactivate**（理由を `history` に追記）。集団全体は保持され、
  `show`/`export_data` で除外理由まで追跡できる。
- **`in_label=True` が既定**: エネルギー処理・`only_minimum` はラベル（=分子）ごとに行う。
  ラベル未設定の System は `""` ラベルの 1 グループとして扱われる。
- **前提処理の自動実行**: `rmsd_limit`→`calc_symm(calc_all=False)`→`calc_bonds`、`write_mol`→`calc_bonds`、
  `get_average`→`calc_distribution`、`GauBox.check_bonding`→`read_vibration`。
- **`Box(box)` と `.plugin.*` は contents を共有（エイリアス）**する（`boxcore.py:39`）。`add(box)` は
  active な System を**参照のまま**追加する。

## 4. Template 機構

実装: `formats.write_input`（`formats.py:205`）→ `replace_key`（`:224`）→ `replace_arg`（`:275`）。

1. テンプレートを `resolve()` し、存在しなければ `ValueError`。System ごとに毎回ファイルを読む。
2. `replace_key`: 文字列全体に対し、下記キーワードが無くなるまでループで置換（置換順は固定）:

| キーワード | 置換値 | 注意 |
|---|---|---|
| `#NAME#` | `c.name` | |
| `#DATA[key]#` | `c.data[key]` | **値が str でないと TypeError**（quirk テストで固定） |
| `#AXYZ#` | `"{:<2} {:>15} {:>15} {:>15}"`（記号, `float_to_str(x/y/z)`） | 精度は Python の repr（丸めなし） |
| `#ATOMS#` | 原子数 | |
| `#CHG#` | `c.charge`（原子形式電荷からの導出を含む） | |
| `#MULT#` | `c.multiplicity`（既定 1） | |
| `#ENRGY#` | `str(c.energy)`（kcal/mol, 未設定なら `"None"`） | 単位・種類は暗黙 |
| `#PATH#` | 書き出し**前**の `c.path` | |
| `#LABEL#` | `c.label` | |

3. `arg: dict` があれば `#KEY#` を `str(value)` で置換（ユーザー定義キーワード）。
4. 出力先: `change_dir(c.path, directory, c.name)` + テンプレートの拡張子。ディレクトリは自動作成。改行は LF 固定。
5. `link=True`（既定）なら `c.path` を新しい入力ファイルへ付け替え → filetype が入力形式に再判定され、
   続く `.run()` / `.submit()` がその入力を実行対象として振り分ける。

その他の性質: 未知の `#XXX#` は黙って残る。置換値にキーワードが含まれると再展開される（自己参照で無限ループの
可能性）。テンプレート内容・パラメータ・展開結果はどこにも記録されない。`TxtBox.parse_keys`
（`txtlib.py:77`）は同じ `replace_key` を使うが、行を連結してしまう（quirk テストで固定）。
PBS ジョブスクリプトも `write_input` で生成でき（`#PBS` を含む `.sh` は `app/pbs/jobscript` と判定）、
`Box.submit()` → `qsub` へ振り分けられる。つまり「計算入力テンプレート」と「実行スクリプトテンプレート」が
同じ機構で扱われている。

## 5. Calculation execution の構造

| 対象 | 実装 | 挙動 |
|---|---|---|
| Gaussian | `gaulib.run/submit`（`gaulib.py:370,382`） | `run`: `subprocess.run([g16, path], cwd)` でブロック。`submit`: `Popen` を投げっぱなし |
| ORCA | `orclib.run/submit`（`orclib.py:85,101`） | 同上。stdout を PIPE で受けてログへ出すのみ（`.out` への書き出しは ACCeL 側では行わない。第2引数の意味は未検証） |
| PBS / OpenPBS | `pbslib.que_submit/que_wait`（`pbslib.py:12,24`） | `qsub` の stdout を `data["jobid"]` に保存。`qstat` を 10 秒間隔でポーリングし、一覧から消えたら完了扱い。タイムアウト・失敗判定なし |
| Slurm | **未実装** | |
| Local 並列 | なし（`submit` の Popen が事実上の並列） | |
| script generation | 専用機構なし（テンプレートで代替） | |
| retry / error handling | `GauBox.check_resub`（`gaulib.py:181`） | ログのエラーシグネチャ分類（link9999 振動判定、SCF 収束失敗、メモリ等）、最終構造から `.gjf` 再生成、旧ファイルを `_` 付きで rename、`data["resubmission"]` |
| コマンド解決 | `Execmd`（`execmd.py`） | キー → 実行ファイルのグローバルマップ |

問題点: `subprocess.run` に `check=True` がないため `CalledProcessError` は発生せず、**失敗が検出されない**。
`submit` はプロセスハンドルを保持しないため完了を追跡できない。実行と出力の対応付けは
「同名ファイルを次のスクリプトで Box に読む」という**ユーザーの運用**に委ねられている。
逆に言えば、v0.3 の典型的な使い方はすでに「既存出力の解析」が中心で、実行は薄い層にとどまっている。

## 6. Parser 構造

- すべて行走査型。メソッド呼び出しごとにファイル全体を再読込。
- 結果を **System に直接書き込み**（`atoms`, `energy`, `charge`, `multiplicity`, `data[...]`）、失敗時は
  `deactivate`。**解析・検証・選別が同じ処理の中で結合**している。
- エネルギーは読み込み時点で kcal/mol に変換して `energy` に上書き。Gaussian は生値（hartree）を
  `data["g16_scf"]` に残すが、ORCA・xTB は生値を残さない。

| filetype | 判定（metric） | 登録された処理 |
|---|---|---|
| `format/xyz` | 拡張子 .xyz (10) | read_atoms |
| `format/mol` | .sdf/.mol/.sd (11) | read_atoms（`data["sdf"]` に生ブロック保存） |
| `format/mae` / `maegz` | .mae (15) / .maegz (20) | read_atoms, read_energy（`r_mmod_Potential_Energy*`, kJ/mol） |
| `app/g16/input` | .gjf/.com/.inp + route 行 (30) | read_atoms（Z-matrix 非対応）, run, submit |
| `app/orca/input` | .com/.inp + 先頭 `!` (40) | run, submit |
| `app/g16/output` | .log/.out/... + "Entering Gaussian System" (50) | read_atoms(archive / standard / input orientation), read_energy(SCF, method 指定可), read_thermal, check_end, check_freq, calc_free_energy |
| `app/xtb/output` | .log/.out/.xtb + バナー (50) | read_energy, read_thermal, check_end |
| `app/pbs/jobscript` | .sh/.qsh/.qsub + `#PBS` (50) | submit, run(submit+wait) |
| `app/orca/output` | .log/.out + バナー (60) | read_atoms（同名 .xyz）, read_energy, check_end |

Gaussian 専用（Box からは `.plugin.gau`）: `read_vibration`, `check_bonding`, `read_nmr`, `read_ecd`,
`read_coupling`, `check_resub`/`get_resub`, `get_irc`, `get_trj`, `is_input/is_output`。
ORCA: `check_optimized`。xTB: `read_free_energy`。Maestro: `read_maegz`, `get_unzip`, `write_mae`（外部 `sdconvert`）。

## 7. 主要アルゴリズムと配置先の提案

| アルゴリズム | 現在地 | 性質 | v2 での配置 |
|---|---|---|---|
| 相対エネルギー | `BoxCore.calc_rel_energy` | 集団依存・**energy を上書き** | 独立 algorithm（純関数）＋ Flow annotation（絶対値は保持） |
| energy filtering | `BoxCore.energy_limit` | 集団依存選別（ΔE ≥ threshold を除外、max_limit） | 独立 algorithm → Flow selection operation |
| Boltzmann 分布 | `BoxCore.calc_distribution` | 集団依存 | 独立 algorithm → Flow analysis operation |
| 最小選択 | `BoxCore.only_minimum` | 集団依存選別（同値は先勝ち） | 独立 algorithm → Flow selection operation |
| Boltzmann 平均 | `BoxCore.get_average` | 集約（新 System 生成） | 独立 algorithm → Flow reduction operation |
| RMSD 重複除去 | `topology.rmsdpruning` | 集団依存・**エネルギー順依存・近似探索**（`redundant_check`） | algorithm を「判定を返す」形に分離 → Flow selection |
| 対称 RMSD 行列 | `Modeler.get_symmetry_matrices` | 構造単体から導出 | 独立 algorithm（Structure の派生量として memo 化可） |
| 結合推定・芳香族 | `Modeler.calc_bonds/aromatize` | 構造単体から導出（閾値依存） | 独立 algorithm（結果は Structure の perception 属性） |
| 立体 | `Modeler.calc_stereo` | 構造単体（「開発中」とログに明記） | 独立 algorithm |
| 原子番号マッピング | `topology.map_numbers`, `Modeler.get_maps` | 参照構造依存・**原子順を in-place 変更** | 独立 algorithm（新 Structure を返す transform） |
| 幾何計測 | `Atoms.get_length/angle/dihedral` | 構造単体 | Structure のメソッド（読み取り専用） |
| 幾何操作 | `Modeler.set_length/mirroring/merge/incorporate` | 構造単体 | 新 Structure を返す transform |
| 状態管理 | `System.deactivate`, `Systems.has_*` | ワークフロー状態 | Flow（選別結果は集合として表現）/ Box 互換層 |

## 8. 技術的負債・確認済みの問題（今回は修正していない）

`quirk` マーカー付きテストで現状を固定したもの（★）と、コード読解で確認したもの。

- ★ `System.duplicate()` が `filetype`・`history`・`distribution` を複製しない（`systems.py:168`）。`Box.duplicate()` も同様。
- ★ `Box(box)` / `.plugin.*` が contents を共有する（エイリアス）。
- ★ `search()` を引数なしで呼ぶと `UnboundLocalError`（`boxcore.py:288`）。
- ★ `energy_limit` などはエネルギー未設定の active System があると `TypeError`。
- ★ `#DATA[key]#` は非 str 値で `TypeError`、`#ENRGY#` は未設定で `"None"`。
- ★ `TxtBox.parse_keys` が行を連結する（`txtlib.py:77`）。
- `write_xyz` / `convert_to_mirror` の centering は重心を入力座標の小数桁数で丸めるため近似的（テストで固定）。
- `run`/`submit` が失敗を検出しない（`check=True` なし）。`submit` はハンドルを保持しない。
- `Atom.show()` が z に y を表示（`atoms.py:174`）。
- `is_mae_format`/`is_maegz_format` の `p.suffix not in (".mae")` はタプルでなく文字列包含判定（`maelib.py:139,146`）。
- `read_xyz` の失敗理由が "orca from xyz file"（`formats.py:18`、コピペ由来）。
- `TxtBox.delete_lines_by_key` はキーワードを**含む**行を残す（名前と逆）。`insert_lines` のスライスは要確認（`txtlib.py:52`）。
- `check_freq` は各振動ブロックの先頭 3 本のみ参照し、ブロックごとに判定を繰り返す。
- `calc_symm(calc_all=False)`（`rmsd_limit` の既定経路）はラベル代表の対称行列を同ラベル全配座へコピーし、
  **原子順が同一であることを暗黙に仮定**する（`boxcore.py:387-395`）。
- `rmsdpruning` は `redundant_check=3` による早期打ち切りを持つ近似探索で、結果がエネルギー順に依存する。
  `include_numeric_isomers` は Box から指定できない。`all_perturbation` の既定値が Box（False）と topology（True）で異なる。
- `energy` の種類（SCF / G / 相対値）と単位が記録されない。`calc_rel_energy` は絶対値を破壊する。
- `multiplicity` の既定値 1、`charge` の原子形式電荷からの導出は、読み込み元に情報がない場合も値が入る。
- 未使用・スタブ: `modeler._cal_sym_rmsd`（`Atoms.data` を参照しており実行不能）, `Modeler.incorporate`,
  `Modeler.get_mapped`, `_superimpose`, `_map_to`, `topology.set_chirality`, `gaulib.read_multiple_xyz`, `Atoms.mw`（誤値を返すとログに明記）。
- `Log`: import 時に DEBUG レベルの stderr ハンドラを登録し、終了時に入力ディレクトリへ `.mcl` を書き出す
  （ユーザーデータのディレクトリへの副作用）。`Log.input_dir` は `System()` 生成のたびに書き換わるグローバル状態。
- import 副作用に依存するグローバルレジストリ（`FileType`, `Selectors`）。plugin の import 順が判定順に影響しうる。
- `Elements` は pickle（`elements.pkl`）。差分レビュー不能で、データの由来はコメントのみ。
- README の表と実装のずれ（`write_xyz(change_path=...)` は実際は `link`、`write_input` の `arg` 未記載）。
- `.gitignore` が `tests/*` を無視していた（今回 `tests/baseline/` のみ例外化）。
- `.pre-commit-config.yaml` の固定バージョンが現行環境で動かない（Python 3.14 で確認）: black 21.11b1 は
  click 8.1 以降で `ImportError: cannot import name '_unicodefun'`、pyupgrade v2.29.1 は `ast.Str` 削除により
  `AttributeError`。check-yaml / end-of-file-fixer / trailing-whitespace / mixed-line-ending は動作する。
  バージョン更新は整形差分を生みうるため今回は変更していない。

---

## 9. Box から v2.0 へ引き継ぐもの

### 9.1 「Box らしい書き心地」を構成している具体的要素

1. **入口が 1 つ**: `Box("*.log")` が glob・ディレクトリ・パス列・System 列・Box を同じように受ける。
2. **動詞＋ドメイン語彙のメソッド名**: `read_*` / `check_*` / `calc_*` / `write_*` / `*_limit` / `only_minimum`。
   workflow 用語（task, node, connect, submit graph）が一切現れない。
3. **全メソッドが自身を返す chain**: 記述順＝処理順で、上から読める。
4. **削除ではなく deactivate**: 除外された構造も理由付きで残り、`show()` / `export_data()` で監査できる。
5. **ラベル＝分子の暗黙グルーピング**（`in_label=True` 既定）: ループを書かずに「分子ごとの配座集合」を扱える。
6. **filetype 自動振り分け**: `read_energy()` が Gaussian / ORCA / xTB / Maestro を区別せず動く。
7. **科学的に妥当な既定値**: 3 kcal/mol、0.01 Å、298.15 K。
8. **前提処理の自動実行**: `rmsd_limit()` が結合・対称を自分で用意する。
9. **テンプレート＝計算定義**: 計算内容はユーザーのテンプレートファイルにあり、ACCeL は `#KEY#` を埋めるだけ。
10. **副次出力を chain の途中に置ける**: `write_xyz()` / `export_data()` / `count()` が流れを止めない。
11. **判断理由のログ**: 除外理由・最小エネルギー等が逐一ログに出る。

v2 の Flow では 1〜10 を維持し、11 は provenance として構造化して保存する形へ発展させるのが自然です。

### 9.2 Flow へそのまま（同じ名前・同じ意味で）引き継げる API

選別・解析・計測で、「入力集合 → 出力集合/注釈」として純粋に定義できるもの:

`energy_limit(threshold, max_limit, in_label)`, `only_minimum(in_label)`, `rmsd_limit(...)`,
`calc_distribution(in_label, temperature)`, `calc_energy(keys, unit)`, `calc_free_energy()`,
`labeling(separator, index_list)`, `set_label`, `set_data`, `calc_length/angle/dihedral`,
`calc_bonds`, `calc_symm`, `count`, `show`, `export_data`, `read_atoms`, `read_energy`, `read_thermal`,
`check_end`, `check_freq`（既存出力の解析として）。

### 9.3 名前は維持できるが内部意味論を変えるべき API

| API | Box での意味 | Flow での意味（案） |
|---|---|---|
| `write_input(template)` | ファイルを書き、`path` を付け替える | Template を用いた**Task の宣言**（遅延）。入力生成は Task の一部 |
| `read_*` | System を上書き | 出力 artifact を解析して**新しい版の Structure / 結果**を作る |
| `calc_rel_energy` | `energy` を上書き | 絶対値を保持したまま相対値を注釈 |
| `energy_limit` / `only_minimum` / `rmsd_limit` | `state` を False にする | 出力集合を返す selection（除外は provenance に残る） |
| `run` / `submit` | 即時実行 / 投げっぱなし | 不要になるか、「この Task を実行対象にする」宣言。実行は `start()`（仮称）時に Plan から |
| `labeling` | ラベルを書き換え | グルーピングキーの定義 |
| `convert_to_mirror` / `modify_length` / `map_numbers` | 座標・原子順を in-place 変更 | 新 Structure を生成する transform |
| `get()` | active の Systems | 評価済み結果へのアクセス |

### 9.4 Box 互換層にのみ残すべき API

`bind`, `contents`/`data` 属性への直接アクセス, `set_state`（手動の再活性化）, `duplicate`（immutable な
Flow では不要）, `Box(box)` のエイリアス挙動, `search`, `copy_files`, `zero_fill`（ファイル名の整形）,
`.plugin.*` アクセサ, `PbsBox.wait`, `submit`（投げっぱなし）, `check_resub`/`get_resub` のファイル rename,
`TxtBox` の行編集, `Dialog`。

---

## 10. Structure

### 10.1 Structure へ移せる System 機能

元素列・座標（`Atoms`）、原子形式電荷、全電荷（明示値と原子からの導出を区別して保持）、多重度（**未指定を
表現できるように**）、結合（入力由来か perception 由来かを区別）、幾何計測（length/angle/dihedral）、
座標変換（mirror / set_length / merge）を「新しい Structure を返す」形で。

### 10.2 Structure（またはそれに付随する結果レコード）へ移すべき Box 内部データ

- `energy` → 種類・単位・由来（どの計算・どの parser）付きの複数エネルギー（例: SCF, ZPE 補正, G 補正）。
- `data` のうち計算結果に当たるもの（`g16_scf`, 熱補正, 振動数, 振動モード, NMR, ECD, 結合定数）→ 由来付き結果。
- `path` → 入力/出力 artifact への参照（パスでなく内容ハッシュ＋保存場所）。
- `name` → 表示名（同一性とは分離）。
- `label` → Structure 自体ではなく、Flow のグルーピング用メタデータとして持つか要決定（§8 論点 Q3）。

### 10.3 Structure へ入れるべきでないもの

`state`/`history`（Flow の選別状態）、`distribution`・相対エネルギー（集団に依存）、ラベル単位の操作、
ファイル I/O、テンプレート展開、計算実行、filetype 判定と振り分け、ログ、アルゴリズムの一時領域
（`atom.cache["isomeric_subs_list"]` 等）、対称行列（派生量として別に memo 化）、`Modeler` のアルゴリズム
（Structure を受け取る関数として外に置く）。

---

## 11. Flow

### 11.1 現行コードのうち Flow operation として再利用できる処理

- `energy_limit`, `only_minimum`, `calc_distribution`, `calc_rel_energy`, `calc_energy`, `get_average` の中核ロジック
  （ラベルごとの min、閾値判定、Boltzmann 重み）→ `(id, energy, group)` を受けて保持集合/値を返す純関数に抽出可能。
- `topology.rmsdpruning` → 「どれを除外するか」を返す関数へ（現状は `deactivate` を直接呼ぶ）。
- `Modeler.calc_bonds/aromatize/get_symmetry_matrices/calc_stereo`, 幾何計測・操作。
- `formats.replace_key/replace_arg` → Structure + Template + params → テキストの純関数（レンダリング）。
- 各 parser → artifact → 結果レコードの関数（System への書き込みと deactivate を外す）。
- `check_end`/`check_freq`/`check_resub` のエラー分類 → validator。

### 11.2 Template を Flow task として利用するために必要な変更

1. **レンダリング・書き出し・付け替えの分離**: 現在の `write_input` は 3 つを同時に行う。
2. **Template オブジェクト**（名称未定）: 元パス、内容（バイト列）とそのハッシュ、使用しているプレースホルダ
   の一覧（テキストから静的に抽出可能）、`arg` パラメータ。
3. **依存の宣言**: テンプレートが使うプレースホルダから、Task が Structure のどの属性に依存するかが決まる
   （例: `#ENRGY#` を使うテンプレートならエネルギーも入力 fingerprint に含める必要がある）。
4. **厳格モード**: 未知の `#X#`、非 str の `#DATA[]#`、`None` 値をエラーにする選択肢（既存 Box の挙動は変えない）。
5. **出力側の束縛**: どのプログラム・どの parser で結果を読むか、期待される出力ファイル。
6. **実行テンプレートとの区別**: 計算入力テンプレートと、ジョブスクリプト（PBS 等）テンプレートを別の役割として扱う。
7. **cache invalidation**: 既定はテンプレート内容のバイト完全一致＋パラメータ＋依存属性の fingerprint。
   コメント・空白のみの変更で再計算するかは論点（Q6）。
8. **provenance**: テンプレート内容（またはハッシュ＋保存コピー）、パラメータ、展開後入力のハッシュ、
   出力から読めるプログラムバージョン（Gaussian archive の `Version=` 等）を記録。

### 11.3 増分再計算を阻害している現在の密結合

1. **System の in-place 変更**（`energy` 上書き、原子順の入れ替え、`path` 付け替え）: 版がなく、前回と今回を比較できない。
2. **選別結果が System 上のフラグ**（`state`）: 選別ノードの出力を集合として保持・比較できない。
3. **`write_input` による path 付け替え**: 同じオブジェクトが「出力構造」から「入力ファイル」に変わり、同一性がパスに縛られる。
4. **エイリアス共有**（`Box(box)`, ビュー, plugin）: 変更の影響範囲が追跡できない。
5. **parse・validate・select の結合**: parser が `deactivate` を直接呼ぶ。
6. **実行と結果の未対応付け**: `submit` はハンドルなし、`run` は失敗未検出、出力の読み込みはユーザー運用。
7. **グローバル状態**: `Log.input_dir`, `Execmd`, レジストリ。
8. **選別が単調でない**: `energy_limit` の基準（最小値）や `rmsd_limit` の比較順は集団全体に依存するため、
   新しい低エネルギー配座 D の追加で、既存の A/B/C の選別結果が変わりうる（除外への変化も起こる）。

### 11.4 Structure 単位の差分追跡に必要と思われる情報

- **内容 fingerprint**: 元素列（原子順込み）、座標（丸め方針要決定）、全電荷、多重度、（ユーザー指定の）結合、
  将来的に周期境界・格子、フラグメント定義、拘束。
- **由来 ID（lineage）**: 元 artifact の内容ハッシュ＋ファイル内インデックス（maegz・IRC・trajectory は 1 ファイル多構造）、
  親 Structure、生成した Task。
- **Task identity**: Template 内容ハッシュ、パラメータ、依存する Structure 属性、プログラム名/バージョン、
  （結果に影響する）実行環境設定。
- **Result identity**: 出力 artifact の内容ハッシュ＋parser の種類/バージョン。
- **グループキー**（label）とエネルギーの種類・単位・由来（選別の入力）。

### 11.5 provenance / cache / invalidation の設計上の論点

§12 の Q1〜Q10 を参照。特に重要なもの:

- 同一性を「内容」で定義するか「由来」で定義するか（同一座標の別由来構造をまとめるか）。
- 座標の数値比較（ファイル由来の 10 進表現 vs 計算で生じた浮動小数点）。
- parser の更新で再解析は必要だが再計算は不要、という区別。
- 失敗した計算の cache（決定的失敗と一時的失敗の区別）。
- cache の保存場所と、絶対パスに依存しない可搬性。
- 選別の非単調性（除外に変わった構造の下流結果は削除せず「今回の出力に含まれない」として保持）。

### 11.6 パラメータ変更時に差分だけ伝播させるための内部境界

```text
User-facing Flow chain（Box 語彙）
    ↓ 各メソッドは汎用 operation を 1 つ追加した新しい Flow 定義を返す（immutable）
Operation 定義
    - 集合演算（selection / grouping / analysis / reduction）: 安価。毎回全体を再評価してよい
    - 構造ごとの写像（Template Task / parse / 構造単体の派生量）: 高価。(task key, structure key) で memo 化
    ↓
依存グラフ（op 列＋入出力）
    ↓
Incremental evaluator
    - 各 op の出力を fingerprint（id 集合＋注釈値）で比較し、変化がなければ下流を止める（early cutoff）
    - 写像 op は「今回の入力 id 集合 − cache 済み key」だけを新規 work とする
    - 選別の結果は上流の結果に依存するため、静的な一括計画ではなく段階的（動的）評価になる
    ↓
Execution plan（今回本当に必要な Task の列）
    ↓
Executor interface（submit → future → artifact）── Local / Parsl(Local, Slurm, PBS)
    ↓
Result store（内容アドレスの artifact ＋ task key 索引、削除しない）
```

要点は **「安価な集合演算は毎回全部やり直し、高価な構造ごとの計算だけを構造単位で再利用する」** という分割です。
これにより、`energy_limit` 3→5 で A,B,C,D,E になれば D,E だけが新規 Task になり、`only_minimum` が引き続き A を
選べば出力 fingerprint が変わらないので下流の高コスト計算は再実行されません。ユーザーはこの区別を意識せず、
Box と同じ chain を書くだけで済みます。

---

## 12. Execution / Parser

### 12.1 Parsl へ委譲できる現行処理

プロセス起動（`run`/`submit`）、PBS への投入と `qstat` ポーリング（Parsl の Provider/Executor が担当）、
Slurm 対応（新規）、並列度・ノード確保・ブロック管理、一時的障害に対するリトライ、完了待ち（future）。

### 12.2 ACCeL 側へ残すべき workflow ロジック

何を実行するかの決定（planner）、入力レンダリング、プログラム/parser の選択、出力の検証（`check_end`,
`check_freq`）、**科学的なエラー分類と再投入方針**（`check_resub`: 構造を変えて再投入するのはドメイン知識）、
cache key と provenance、選別の意味論、グルーピング、作業ディレクトリと artifact の命名規則、実行コマンドの解決（`Execmd` 相当の設定）。

### 12.3 cclib へ移行可能な parser 処理（要パリティ検証）

Gaussian: 座標（archive/orientation に相当する `atomcoords` の最終構造）、原子番号、電荷・多重度、SCF エネルギー、
正常終了判定、振動数（`vibfreqs`）、振動モード（`vibdisps`）、熱化学量、励起状態（`etenergies`, `etoscs`, `etrotats`）、
最適化軌跡（`get_trj` 相当）。ORCA: 最終エネルギー、座標（同名 .xyz 依存を解消できる）、正常終了。

注意: cclib の**単位規約**（エネルギー単位はバージョンにより異なり得る）と、熱化学量が「補正値」ではなく
「電子エネルギー＋熱補正」の合計で提供される点など、ACCeL の `g16_corr_to_gibbs` と定義が異なる可能性がある。
**同じ出力に対して既存 parser と cclib の結果を突き合わせるパリティテストを先に作ること**が移行の条件です。
xTB・NMR テンソル・スピン-スピン結合の cclib 対応状況は未確認（要調査）。

### 12.4 ACCeL 独自 parser として残すべき処理

`check_resub`（エラー分類＋入力再生成）、`check_bonding`（振動変位順位のヒューリスティック）、
method 指定付き SCF 抽出（`SCF Done:  E(method)`）、入力ファイルからのジョブ数推定（`--link1--`, opt+freq）、
Maestro `.mae/.maegz`、sdf/mol（生ブロック保持を含む）、Gaussian 入力 `.gjf` の読み込み、xyz、
filetype 判定（入力ファイル・ジョブスクリプトの判定は cclib の範囲外）、IRC の方向・反応座標、
cclib が提供しない場合の結合定数の成分選択・xTB の熱化学内訳。

### 12.5 ASE が有効と思われる用途

ファイル形式の相互運用 fallback（extxyz, cif, POSCAR, pdb 等）、特に**結晶・周期系**（Structure の一般化で必要になる
格子・周期境界）、`ase.Atoms` との相互変換による外部ツール連携、大規模系の近傍探索。単位・原子順・周期情報の
変換規則は明示的にテストすること。

---

## 13. 移行ロードマップ（提案）

最初に設計すべきは **(a) Structure の同一性/fingerprint の仕様** と **(b) Operation 契約（安価な集合演算と
高価な構造ごとの写像の区別、出力 fingerprint による early cutoff）** の 2 点です。Box と同じ chain を保ったまま
増分再計算を成立させる鍵はこの 2 つであり、`Structure` のクラス設計も Flow のメソッド設計もここから導かれます。

| Phase | 内容 | 完了条件 |
|---|---|---|
| 0 | baseline tests / CI / CLAUDE.md / 本監査（**今回**） | 既存挙動が CI で固定されている |
| 1 | 同一性・fingerprint・Operation 契約の設計文書と、§14 の論点への決定 | 決定記録（ADR）がある。実装なし |
| 2 | `Structure`（immutable 値オブジェクト）と `System ⇄ Structure` 変換アダプタ | Box 無変更、往復変換テスト |
| 3 | アルゴリズムの純関数化（energy/Boltzmann/only_minimum/RMSD 判定を「判定を返す」形に）。Box はそれを呼ぶだけに | baseline tests が不変のまま通る |
| 4 | parser 層: 「artifact → 結果レコード」関数へ分離し、cclib backend を追加して既存 parser とのパリティテスト | 同一出力で値が一致（差分は文書化） |
| 5 | Template オブジェクト（純粋レンダリング、プレースホルダ抽出、fingerprint）。`write_input` は内部でこれを使う | 出力バイト一致 |
| 6 | Flow 定義モデル（immutable chain、Box 語彙）と、**計算を伴わない解析専用 Flow** のインメモリ評価器 | 既存出力だけで Box と同じ結果 |
| 7 | 増分評価エンジン: result store、構造単位 memo、early cutoff、段階的評価。**決定的な fake executor** でテスト | 3→5 kcal/mol、最小構造不変/変化のシナリオが自動テストで検証される |
| 8 | Executor インタフェースと Local 実装 → Parsl backend（Local/Slurm/PBS） | fake と同じシナリオが Local で通る |
| 9 | Box 互換統合（Box 内部で Structure/algorithm を利用）と Flow UX の磨き込み | baseline tests 不変 |

例示された順序からの主な変更点: 同一性仕様を Structure 実装より前に置く（Phase 1）、parser 層の分離を
Template/Flow より前に置く（解析専用 Flow が parser の純関数化に依存するため）、Parsl は fake executor で
エンジンを完成させた後に接続する（エンジンの正しさをスケジューラから切り離して検証するため）。

---

## 14. 決定が必要な論点（推測で確定していないもの）

- **Q1 Structure の同一性**: 内容（元素・座標・電荷・多重度）で定義するか、由来（元ファイル＋インデックス＋親）で
  定義するか、両方を持つか。同一座標だが別由来の構造を 1 つにまとめてよいか。
- **Q2 座標の比較精度**: fingerprint の座標丸め（例: 1e-6 Å）を使うか、ファイル由来の 10 進文字列を正とするか。
  丸め境界での誤判定をどう扱うか。
- **Q3 label の位置付け**: Structure の属性か、Flow のグルーピング定義か。v0.3 はファイル名から導出。
- **Q4 原子順**: 原子順だけが異なる同一分子を同一とみなすか（計算入力としては別物、比較では同一）。
- **Q5 エネルギーの種類**: `energy_limit` 等がどのエネルギー（SCF / G / 相対）を使うかを Flow でどう指定・記録するか。
  v0.3 は「最後に書かれた値」。
- **Q6 Template の変更検出**: バイト完全一致か、空白・コメント差を無視するか。
- **Q7 parser 更新の扱い**: parser バージョン変更時に再解析のみを行う仕組みを持つか。
- **Q8 失敗の cache**: SCF 未収束などの決定的失敗を cache するか、再試行するか。`check_resub` 相当の自動再投入を Flow の既定にするか。
- **Q9 パラメータ変更の UX**: (a) スクリプトの値を書き換えて再実行するだけで永続 cache が効く方式（Box の使い方に最も近い）、
  (b) 名前付きステップとパラメータ差し替え API、のどちらを主とするか。本監査は (a) を主、(b) を補助とすることを推奨。
- **Q10 Flow における実行の明示**: `write_input(template)` の後に `run()` を要求するか、`write_input` → `read_*` の並びで
  Task を暗黙に成立させるか。
- **Q11 既存 quirk の扱い**: §8 の★項目（`duplicate` の欠落、`search()` 等）を v0.x で修正するか、v2 の Box 互換層でも保持するか。
- **Q12 multiplicity 既定値 1 / 原子形式電荷からの全電荷導出** を Structure でも既定とするか、「未指定」を必須にするか。
