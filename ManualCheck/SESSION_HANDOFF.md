# 引き継ぎメモ: Pickup_reference_track_for_each_track / Check_track_consistency

このドキュメントは、NINJA E71a実験のManualCheckパイプライン（`C:\Users\kasumi\source\repos\ManualCheck\`）のうち、
`Pickup_reference_track_for_each_track` と `Check_track_consistency` の2ツールについて、
前セッションまでの作業内容をまとめたものです。新しいセッションでこの続きから作業する際は、
このファイルの内容を最初に読み込んでください。

## プロジェクト概要

- NINJA E71a実験の原子核乾板飛跡データを、手動チェック（マニュアルチェック）用にリストアップ・検証するツール群
- `Pickup_reference_track_for_each_track.cpp`: TrackList(YAML)を読み、各飛跡（実在/外挿）に対してReference（比較用飛跡）を1つ選び出力する
- `Check_track_consistency.cpp`: 上記の出力が、Area0(merge)/Area1~6(小エリア)/fvxxの3種のvxxファイル間で矛盾なく追跡できるかを検証する独立ツール（今セッションで新規作成）
- データはvxxフォーマット（`VxxReader.lib`、ソースは `E:\udd\hayakawa\VxxReader\` にあり、疑問があればここのSample.cpp/VxxReader.cppを直接参照するとよい）
- ECC4のデータは Area0とArea1が `K:\NINJA\E71a\ECC4\...`、Area2~6が `I:\NINJA\E71a\ECC4\...` に分かれて保存されている（このルールはコード中の複数箇所に埋め込まれている）

## 現在の状態

両ツールとも**ビルド成功・想定通りに動作することを確認済み**（クラッシュなし）。

## `Pickup_reference_track_for_each_track.cpp` — 今セッションでの変更点

### 出力形式
- 実在/外挿/Referenceの出力は、以前追加した `#RawID`/`#Track`/`#Zone` コメント行をすべて削除し、元の素の形式に戻した:
  - 実在: `Plate / RawID / Surface / SaveTo`
  - 外挿: `Plate / Track / Zone / Surface / SaveTo`
  - Reference: `Plate / RawID / Zone / Surface / SaveTo`
- 小エリア側の追加情報（`RawID_small`, `Track_small`, `Zone_small`）は、`out_dir_path\crosscheck.yaml` という**単一ファイル**（実行全体で1つ、`PLxxx.yaml`と同じディレクトリ）にまとめて出力するようにした。
- `refSelectMode` 引数（`argv[5]`, 省略可）: 0=ランダム選択(既定), 1=VPHの和(`m[0].ph%10000+m[1].ph%10000`)が最大のものを選択。

### バグ修正
1. 実在トラックのReference探索が、マッピング後の小エリア(`Area{itr->Zone}`)ではなく`Area0`のmergedファイルを検索していた座標系不一致バグを修正（`file_in_ref_area`を新たに構築して渡すように変更）。
2. `argc < 3` のチェックが不十分で、`argv[4]`(output_log)未指定時に `nullptr` から `std::string` を構築してクラッシュしていたのを `argc < 5` に修正。
3. 出力先ディレクトリが既存の場合に実際には停止しない `exit;`（`exit()` の呼び出し忘れ、セミコロンで関数参照だけして終わっていた）を `exit(1);` に修正。
4. 一時的に追加していた「Reference basetrackの`.zone`フィールドが実在/外挿のzoneと一致するか」のチェックは、参照している`base_track_t.zone`フィールドの意味が不明瞭（Area番号とは限らない）で誤検知だったため削除済み。

### クラッシュ対策（重要）
- `Search_basetrack` / `Search_nearest_basetrack` で小エリアbvxxファイル（`Area{1~6}\PL{pl}\b{pl}.sel.cor.vxx`）を読み込む際、特定のファイルで**再現性を持ってクラッシュする**問題が発生した（Windows例外コード `0xC0000409` = スタックバッファのセキュリティチェック違反、`ucrtbase.dll`内、`VxxReader.lib`内部の`PMicroTrackBlock`/`BlockManager`相当の処理と思われる）。
  - 一部のケース（PL016）は入力TrackListデータ自体の誤りが原因と判明（ユーザー確認済み）。
  - しかし別のケース（PL091, RawID=4457418）は入力データが正しいにもかかわらず再現し、Area1~6の総当たり探索中に特定のArea読み込みでクラッシュした。**ライブラリ側の一般的な信頼性問題**と判断。
- 対策として、危険な読み込み処理を**子プロセスに分離**した:
  - `main()` の冒頭に `--readarea <file> <pl> <area> <filtered> <xmin> <xmax> <ymin> <ymax> <outfile>` というワーカーモードを追加。自分自身をこのモードで再起動すると、1回の`ReadAll`呼び出しだけを行いバイナリファイルに結果を書き出して終了する。
  - `SafeReadAreaBvxx()` 関数が `CreateProcess` でこのワーカーを起動し、`WaitForSingleObject`で終了を待ち、終了コードが0以外（＝クラッシュ）なら `[WARN]` ログを出して「そのAreaでは見つからなかった」扱いにして処理を継続する。正常終了なら結果ファイルを読んで通常通り処理する。
  - `Search_basetrack`/`Search_nearest_basetrack`はこの`SafeReadAreaBvxx`経由で読むように書き換え済み。
  - **`vxx::opt::a`（CutAreaフィルタ）は使わず、フィルタなしの全件読み込みを維持**している（過去にフィルタを有効化したところ「見つからないことがある」との理由でユーザーが手動でフィルタなしに戻した経緯があるため、その挙動を尊重）。
- Area1~6のどこにも一致するbasetrackが見つからなかったケース（`Basetrack ... was not found in small area.`）は `notFoundInSmallAreaLog`（`SmallAreaNotFoundEntry`のvector）に記録し、プログラム終了時にまとめて `PL, event, rawid(merged), ax, ay, x, y` を出力するようにした。

## `Check_track_consistency.cpp` — 今セッションで新規作成

- 目的: `Pickup_reference_track_for_each_track.cpp`が出力したPLxxx.yamlを読み、実在/Referenceの各飛跡が下記3ファイル間で矛盾なく追跡できるかを独立に検証する。
  - `Area{zone}\PL{pl}\b{pl}.sel.cor.vxx`（小エリア側 basetrack）
  - `Area{zone}\PL{pl}\f{pl}2_thick_0.vxx`（小エリア側 fvxx, surface=2のみ）
  - `Area0\PL{pl}\b{pl}.sel.cor.vxx`（mergeエリア側 basetrack）
- VSプロジェクト・ソリューション(`ManualCheck.sln`)への登録済み。ビルド成功確認済み。
- **入力ファイルはPickup側の出力形式変更に依存しない**設計にした: `PLxxx.yaml`（Pickup側の素の出力）だけを読み、小エリア側の情報（zone/rawid/Track）は**このツール自身がArea0とArea1~6を独立に読み比べて導出**する（`Pickup...`の`crosscheck.yaml`は読まない・依存しない）。
  - 実在トラック: `VerifyReal()` — Area0のbasetrack(`RawID`で直接検索)のmicrotrackペア(m[0]/m[1].rawid)と一致するbasetrackをArea1~6から探す（`GuessArea`で位置から所属Areaを推定し、まずそこを試し、だめなら他のAreaも順に試す）。
  - Referenceトラック: `VerifyReference()` — YAMLに直接ある`RawID`/`Zone`（既に小エリア側）を使い、fvxx確認とArea0側の逆引き（microtrackペア一致するbasetrackをArea0の全件から探す、PLごとにキャッシュ）を行う。
- **重要な技術的発見**: `FvxxReader::ReadAll`/`Begin`の第2引数は**PL番号ではなく"POS"**であり、`POS = PL番号*10 + face`（face=1or2, fvxxファイル名`f{PL}{face}_thick_0.vxx`のfaceに対応）で計算する必要がある。これはVxxReaderライブラリのソース(`E:\udd\hayakawa\VxxReader\Project\Sample\Sample.cpp`のコメント、および`VxxReader.cpp`の`MakePos`関数)から判明した。当初この計算をしておらず、fvxx検証が全件失敗していた。
- fvxxの一致判定は最終的に「rawidのみで判定」「surface2ファイルのみ確認」「m[0]/m[1]のどちらか片方が見つかればOK（両方は要求しない）」という形に単純化されている。

## ファイルエンコーディングに関する注意（重要）

- `Check_track_consistency.cpp` は**UTF-8 BOM付き**で保存する必要がある（`\xEF\xBB\xBF`で始まる）。BOMが無いとMSVCがコードページ932(Shift-JIS)として誤解釈し、日本語コメントが原因でコンパイルエラーの連鎖が起きる。編集後は必ずBOMが保持されているか確認すること。
- `Pickup_reference_track_for_each_track.cpp` はBOM無しで問題なくビルドできている（日本語コメントの配置が偶然問題を起こさない位置にあるため）。新しくファイルを作る場合は基本的にBOM付きにしておくのが安全。

## 未対応・今後の検討事項

- `Search_ref_btrk`（Reference候補探索）も同様の`br.ReadAll`パターンで小エリアbvxxを読んでいるが、まだプロセス分離していない。今後クラッシュが報告されたら同様の対策を検討する。
- `Check_track_consistency.exe`は現状、危険な読み込み処理（Area1~6の全件スキャン等）をプロセス分離していない。同様のクラッシュが起きる可能性があるので、報告があれば`Pickup...`と同じ手法を適用する。
- クロスチェック用ファイル`crosscheck.yaml`は現状誰も自動検証していない（人間が目視確認する想定）。

## 参考: 使い方（現状のCLI引数）

```
Pickup_reference_track_for_each_track.exe #ECC [in-YAML-dir] [output_dir] [output_log_path] [refSelectMode(0=random,1=maxVPH, optional)]
Check_track_consistency.exe #ECC [Pickup...の出力ディレクトリ(PLxxx.yamlがある場所)] [output-log-file-path]
```
