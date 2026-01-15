%[text] ## 解析用のデータを作成
%[text]  解析途中で基準を変更して解析したい場合がある。
%[text] オリジナルのデータをそのまま解析していたら，最初からデータを作成しなければならないので，解析データとオリジナルデータは分けること。

% MMN_raw から年度ごとのデータを読み込み
if ~exist('inputDir_mmn', 'var') || isempty(inputDir_mmn)
    thisFileDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisFileDir)); % .../脳計測科学
    inputDir_mmn = fullfile(projectRoot, 'MMN_raw');
end

% 利用可能な年度を検出
mmn_files = dir(fullfile(inputDir_mmn, '*data_mmn.mat'));
if isempty(mmn_files)
    error('MMN_raw フォルダにデータが見つかりません: %s', inputDir_mmn);
end

available_years = [];
for i = 1:length(mmn_files)
    tokens = regexp(mmn_files(i).name, '^(\d{4})data_mmn\.mat$', 'tokens');
    if ~isempty(tokens)
        available_years(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
    end
end
if isempty(available_years)
    error('MMN_raw フォルダに有効な年度データが見つかりません');
end
available_years = sort(available_years);
fprintf('利用可能な年度: %s\n', num2str(available_years));

% 閾値（未設定なら80）
if ~exist('Amplifired_threshold', 'var') || isempty(Amplifired_threshold)
    Amplifired_threshold = 80;
end

% 出力先の準備
if ~exist('outputDir_filtered','var') || isempty(outputDir_filtered)
    outputDir_filtered = fullfile(projectRoot, 'MMN_filtered');
end
if ~exist(outputDir_filtered, 'dir')
    mkdir(outputDir_filtered);
end

% 全年度を順次処理
for target_year = available_years
    input_file = fullfile(inputDir_mmn, sprintf('%ddata_mmn.mat', target_year));
    fprintf('\n=== %d 年度の処理開始 ===\n', target_year);
    fprintf('データ読み込み中: %s\n', input_file);
    S = load(input_file);
    if ~isfield(S, 'data0') || ~isfield(S, 'trg')
        warning('data0 または trg が存在しません: %s (スキップ)', input_file);
        continue;
    end
    data0 = S.data0; trg = S.trg;
    fprintf('年度 %d のデータサイズ: %s\n', target_year, mat2str(size(data0)));
    data = data0;
    sz = size(data);

    % artifact detection（ベクトル化）
    mx = max(data, [], 2, 'omitnan');
    mn = min(data, [], 2, 'omitnan');
    mask = (mx - mn) > Amplifired_threshold;
    mask_expanded = repmat(mask, [1 sz(2) 1 1 1]);
    data_before_reject = data; %#ok<NASGU>
    data(mask_expanded) = NaN;

    % 保存（閾値入りファイル名）
    outFile = sprintf('%ddata_filtered_th%d.mat', target_year, Amplifired_threshold);
    save(fullfile(outputDir_filtered, outFile), 'data', 'trg', '-v7.3');
    fprintf('保存完了: %s\n', fullfile(outputDir_filtered, outFile));
end
fprintf('\n=== 全年度の処理が完了しました ===\n');
% 今後dataを解析していくこととします。
%%
%[text] ## artifact detection（別の方法）
%[text] 上記とどちらか一方だけの解析でよいです。
% data0->(chl,time,epoch,freq,prt); (31,800,1000,3,3); 
Amplifired_threshold = 80;

sz = size(data);

% 第2次元方向(時間軸など)で最大・最小を求める（NaNは無視）
% epochの最大の部分と最小の部分の値を求めていル。2次元方向で最大の部分
mx = max(data, [], 2, 'omitnan');

mn = min(data, [], 2, 'omitnan');

% 差分が100を超える箇所を判定するマスクを生成

mask = (mx - mn) > Amplifired_threshold; 
% maskのサイズは [sz(1),1,sz(3),sz(4),sz(5)]

% マスクを第2次元方向へ拡張（ブロードキャスト）
mask_expanded = repmat(mask, [1 sz(2) 1 1 1]);

% 元データを退避（除去エポック可視化用）
data_before_reject = data;

% マスクがtrueの箇所をNaNへ置き換える
data(mask_expanded) = NaN;

% %% 除去されたエポックの波形可視化
% % チャンネル方向いずれかで閾値超過したエポックを抽出
% rejected_epoch = squeeze(any(mask, 1)); % (epoch,freq,prt)

% time_idx = 1:size(data, 2);

% % maskから直接超過チャンネルを抽出
% % mask: (ch, 1, epoch, freq, prt)
% % 各チャンネルで少なくとも1つのエポック/条件/参加者で超過があるチャネルを探す
% mask_any_ch = any(mask(:), 'all');  % どのチャンネルでも超過があるか確認
% if mask_any_ch
%     % 各チャネルについて超過があるかをチェック
%     exceeded_per_ch = zeros(size(data, 1), 1);
%     for ch = 1:size(data, 1)
%         exceeded_per_ch(ch) = any(mask(ch, :, :, :, :), 'all');
%     end
%     idx = find(exceeded_per_ch, 1, 'first');
% else
%     idx = 1;  % 超過がなければ1ch
% end
% plot_ch_rej = idx;      % 表示チャンネル

% % デバッグ表示：maskの True 要素をチャンネルごとに表示
% fprintf('\n=== Artifact Detection Debug Info ===\n');
% fprintf('超過閾値: %.1f µV\n', Amplifired_threshold);

% % maskの True 要素のインデックスを取得
% [ch_idx, ~, epc_idx, freq_idx, prt_idx] = ind2sub(size(mask), find(mask));
% n_true = length(ch_idx);

% if n_true > 0
%     fprintf('maskで True となった要素数: %d\n', n_true);
    
%     % チャンネルごとに集約
%     unique_ch = unique(ch_idx);
%     fprintf('\n超過したチャンネル一覧：\n');
%     for uch = unique_ch'
%         count_ch = sum(ch_idx == uch);
%         fprintf('  Ch %d: %d件の超過\n', uch, count_ch);
%     end
    
%     fprintf('\n詳細（最初の30件）：\n');
%     fprintf('(ch, epoch, freq, prt):\n');
%     for i = 1:min(n_true, 30)
%         fprintf('  [%d, %d, %d, %d]\n', ch_idx(i), epc_idx(i), freq_idx(i), prt_idx(i));
%     end
%     if n_true > 30
%         fprintf('  ... (%d件省略)\n', n_true - 30);
%     end
% else
%     fprintf('maskで True となった要素: なし\n');
% end

% % maskで True となっているチャンネルのみをプロット
% fprintf('\nプロット対象チャンネル: %s\n', num2str(unique_ch'));
% max_plot_each = 5;     % チャンネル×dev×prtごとに最大何エポックまで描画するか

% for plot_ch = unique_ch'
%     for prt = 1:size(mask, 5)
%         for dev = 1:size(mask, 4)
%             % このチャンネル・周波数・参加者でmask=Trueとなっているエポック
%             rej_list = find(squeeze(mask(plot_ch, 1, :, dev, prt)));
%             if isempty(rej_list)
%                 continue;
%             end
%             n_show = min(numel(rej_list), max_plot_each);
%             figure; tiledlayout(n_show, 1);
%             for k = 1:n_show
%                 epc = rej_list(k);
%                 nexttile;
%                 plot(time_idx, squeeze(data_before_reject(plot_ch, :, epc, dev, prt)), 'r');
%                 title(sprintf('Rejected epoch %d (dev=%d, prt=%d, ch=%d)', epc, dev, prt, plot_ch));
%                 xlabel('Sample index'); ylabel('Amplitude (\muV)'); grid on;
%             end
%         end
%     end
% end

% %% アーティファクト除去後の平均波形プロット
% % data は (ch,time,epoch,freq,prt)。NaN を含むので mean(...,'omitnan') を使用。
% time_idx = 1:size(data, 2); % サンプルインデックス（実時間が必要ならサンプリング周波数に応じて軸を変換）

% % まず epoch 方向で平均（NaN無視） -> (ch,time,freq,prt)
% avg_ep = squeeze(mean(data, 3, 'omitnan'));

% % 参加者方向でも平均（NaN無視） -> (ch,time,freq)
% grand_avg = squeeze(mean(avg_ep, 4, 'omitnan'));

% % プロットするチャンネルを選択（例: 1ch）。別の電極を見たい場合は適宜変更。
% plot_ch = 14; %fz 5 cz 14

% figure; hold on;
% clr = lines(size(grand_avg, 3));
% for dev = 1:size(grand_avg, 3)
%     plot(time_idx, squeeze(grand_avg(plot_ch, :, dev)), 'Color', clr(dev, :), 'LineWidth', 1.5);
% end

% legend({'1010 Hz','1030 Hz','1050 Hz'}, 'Location', 'best');
% xlabel('Sample index');
% ylabel('Amplitude (\muV)');
% title(sprintf('Artifact-removed grand average (Ch %d)', plot_ch));
% ylim([-2 2])
% grid on;

%% アーティファクト除去済みデータの保存
% 出力先フォルダ設定（ユーザーが `outputDir_filtered` を事前に定義していなければ既定値を使用）
if ~exist('outputDir_filtered','var') || isempty(outputDir_filtered)
    thisFileDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisFileDir)); % .../脳計測科学
    outputDir_filtered = fullfile(projectRoot, 'MMN_filtered');
end
if ~exist(outputDir_filtered, 'dir')
    mkdir(outputDir_filtered);
end

% 年度に応じた出力ファイル名を自動決定（folder パスから 20xx を抽出）
if exist('folder', 'var')
    yrToken = regexp(folder, '(20\d{2})', 'match', 'once');
else
    yrToken = '';
end
if isempty(yrToken)
    outFile = sprintf('data_filtered_th%d.mat', Amplifired_threshold);
else
    outFile = sprintf('%sdata_filtered_th%d.mat', yrToken, Amplifired_threshold);
end

% アーティファクト除去済みの data と trg を保存
save(fullfile(outputDir_filtered, outFile), 'data', 'trg', '-v7.3');
fprintf('\nアーティファクト除去済みデータを保存しました: %s\n', fullfile(outputDir_filtered, outFile));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":9.2}
%---
