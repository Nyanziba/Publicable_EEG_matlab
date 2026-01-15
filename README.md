Public EEG Script
==================

概要
----
授業・自習向けに公開するERP/MMN解析用のスクリプト集です。自分で改変したプログラムは拡張子が .m のファイルです（.mlx は元のライブスクリプトのまま）。MATLAB での実行を想定しています。

リポジトリ構成
----------------
- README.md : この説明書
- practice2025_lecture/ : 授業で配布するスクリプトとデータ
	- MMNanalysis001.m 〜 MMNanalysis008_MIXanova_ans.m : 各回の分析ステップをコード化したファイル（改変済み）
	- *live.mlx : 上記と同内容のライブスクリプト版（未改変）
	- EEG_locs_new.mat : 電極位置情報
	- makeSubplotsClickable.m : サブプロット選択用の補助スクリプト

前提環境
--------
- MATLAB R2022b 以降を推奨（Signal Processing Toolbox, Statistics and Machine Learning Toolbox を使用）
- OS: macOS/Windows いずれも可

使い方（基本の流れ）
--------------------
1. practice2025_lecture/ を MATLAB のカレントフォルダに設定。
2. 目的の .m ファイルを開き、コメントに沿って上から順に実行。
	 - 例: MMNanalysis005.m はフィルタ適用後の波形確認、MMNanalysis007_repANOVA.m は反復測定 ANOVA。
3. 解析途中で生成される図や結果はスクリプト内で指定したパスに保存されます。必要に応じて出力先を変更してください。

各スクリプトの役割（簡易メモ）
------------------------------
- MMNanalysis001〜004: エポック抽出、条件別平均、基本的な可視化。
- MMNanalysis005/005_filtered: フィルタ処理と参加者別の波形・頭皮分布プロット。
- MMNanalysis006: 年度別データからグランドアベレージを生成。
- MMNanalysis007_repANOVA: 単一要因の反復測定 ANOVA の雛形。
- MMNanalysis008_MIXanova_ans: 混合要因 ANOVA 例と事後解析。
- makeSubplotsClickable: サブプロットをクリックで拡大表示する補助。

注意点
------
- .m 版のみコード修正済みです。ライブスクリプト版 (.mlx) は参照用に残しています。
- データファイル (.mat) は実験データを含むため取り扱いに注意してください。
- スクリプトの冒頭にあるパラメータ（フィルタ設定、時間窓、参加者 ID など）は各自のデータに合わせて調整してください。

質問・改善
----------
改善案や不具合があれば、該当ファイル名と事象を添えて連絡してください。
