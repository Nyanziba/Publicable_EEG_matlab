%[text] # フィルタリング済みデータから個人波形をプロット・保存
%[text] MMN_filtered フォルダから各年度のデータを読み込み、参加者ごとにプロットして保存します。
inputDir_filtered = '/Users/nyanziba/school/脳計測科学/MMN_filtered'; % フィルタリング済みデータのフォルダ
outputDir_individual = '/Users/nyanziba/school/脳計測科学/MMN_individual_plots'; % 個人波形の保存先フォルダ
% 入力・出力フォルダの設定
if ~exist('inputDir_filtered', 'var') || isempty(inputDir_filtered)
    thisFileDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisFileDir));
    inputDir_filtered = fullfile(projectRoot, 'MMN_filtered');
end

if ~exist('outputDir_individual', 'var') || isempty(outputDir_individual)
    outputDir_individual = fullfile(projectRoot, 'MMN_individual_plots');
end
if ~exist(outputDir_individual, 'dir')
    mkdir(outputDir_individual);
end

% 利用可能な年度を検出
filtered_files = dir(fullfile(inputDir_filtered, '*data_filtered*.mat'));
if isempty(filtered_files)
    error('MMN_filtered フォルダにデータが見つかりません: %s', inputDir_filtered);
end

available_years = [];
for i = 1:length(filtered_files)
    tokens = regexp(filtered_files(i).name, '^(\d{4})data_filtered', 'tokens');
    if ~isempty(tokens)
        available_years(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
    end
end
if isempty(available_years)
    error('有効な年度データが見つかりません');
end
available_years = unique(sort(available_years));
fprintf('利用可能な年度: %s\n', num2str(available_years));

% 全年度を処理
for target_year = available_years
    fprintf('\n=== %d 年度の個人波形作成 ===\n', target_year);
    
    % フィルタリング済みデータを読み込み
    filtered_file = fullfile(inputDir_filtered, sprintf('%ddata_filtered_th*.mat', target_year));
    files = dir(filtered_file);
    if isempty(files)
        warning('年度 %d のフィルタリング済みデータが見つかりません (スキップ)', target_year);
        continue;
    end
    load(fullfile(inputDir_filtered, files(1).name));
    fprintf('読み込み: %s\n', files(1).name);
    fprintf('データサイズ: %s\n', mat2str(size(data)));
    
    % 参加者ID（Dataset フォルダから取得）
    dataset_folder = sprintf('/Users/nyanziba/school/脳計測科学/MATLAB_EEG/Dataset/%dAutumn_mmn_set', target_year);
    if exist(dataset_folder, 'dir')
        mat_files = dir(fullfile(dataset_folder, '*.mat'));
        names = {mat_files.name};
        prefixes = cellfun(@(s) s(1:8), names, 'UniformOutput', false);
        fn0 = unique(prefixes);
    else
        % フォルダがない場合は参加者番号で命名
        fn0 = cell(1, size(data, 5));
        for p = 1:size(data, 5)
            fn0{p} = sprintf('%d_prt%02d', target_year, p);
        end
    end
    
    %% plot ERP data(individual)
    fs = 1000; % サンプリング周波数
    N = size(data, 2); % データ長
    
    % ----------trigger
    marker = zeros(2, 2);
    marker(:, 1) = 0.0;
    marker(1, 2) = -30;
    marker(2, 2) = 30;
    % ----------trigger
    
    x = linspace(-0.1, (N-100)/fs, N);
    
    load('EEG_locs_new.mat'); % channel位置を取得する。
    
    sz = get(0, 'ScreenSize'); % スクリーンサイズを取得する。
    
    % data->(ch,time,epoch,freq,prt); (31,800,1000,3,X);
    
    for prt = 1:size(data, 5)
        figure('Position', [sz(1) sz(2) sz(3) sz(4)]) % 画面いっぱいのfigureを作成する
        for chl = 1:size(sloc, 1)
            subplot('position', [sloc(chl, 1) sloc(chl, 2) 0.07 0.07])
            hold on
            % ---------------------------- y1-y6をMMNanalysis004.mlxからコピペ
            %findなんてしないで、論理演算子で
            y1 = mean(data(chl, :, trg(:, 1, prt) == 1, 1, prt), 3, "omitnan"); % 1010Hz % standard
            y2 = mean(data(chl, :, trg(:, 1, prt) == 2, 1, prt), 3, "omitnan"); % 1010Hz % deviant
            y3 = mean(data(chl, :, trg(:, 2, prt) == 1, 2, prt), 3, "omitnan"); % 1030Hz % standard
            y4 = mean(data(chl, :, trg(:, 2, prt) == 2, 2, prt), 3, "omitnan"); % 1030Hz % deviant
            y5 = mean(data(chl, :, trg(:, 3, prt) == 1, 3, prt), 3, "omitnan"); % 1050Hz % standard
            y6 = mean(data(chl, :, trg(:, 3, prt) == 2, 3, prt), 3, "omitnan"); % 1050Hz % deviant
            % ----------------------------
            plot(x, y1, 'Color', [255 75 0]/255, 'LineWidth', 2, 'LineStyle', '--'); % 1000Hz % standard
            plot(x, y2, 'Color', [255 75 0]/255, 'LineWidth', 2); % 1010Hz % deviant
            plot(x, y3, 'Color', [0 90 255]/255, 'LineWidth', 2, 'LineStyle', ':'); % 1000Hz % standard
            plot(x, y4, 'Color', [0 90 255]/255, 'LineWidth', 2); % 1030Hz % deviant
            plot(x, y5, 'Color', [3 175 122]/255, 'LineWidth', 2, 'LineStyle', '-.'); % 1000Hz % standard
            plot(x, y6, 'Color', [3 175 122]/255, 'LineWidth', 2); % 1050Hz % deviant
            axis('ij')
            ymx = 10;
            xlim([-0.1 0.7])
            ylim([-ymx ymx])
            title([ ch_name{chl} '(' num2str(chl) ')' ]) % titleを付けます。
            plot(marker(:, 1), marker(:, 2), 'k:', 'LineWidth', 2); % trigger
            set(gca, 'box', 'on');
            set(gca, 'linewidth', 1)
        end
        legend('1000Hz', '1010HzDev', '1000Hz', '1030HzDev', '1000Hz', '1050HzDev')
        makeSubplotsClickable();
        
        % 年度付きファイル名で保存
        outname = sprintf('%d_%s', target_year, fn0{prt});
        saveas(gcf, fullfile(outputDir_individual, outname), 'fig');
        fprintf('保存: %s.fig\n', outname);
        close(gcf);
    end
    fprintf('年度 %d の処理完了: %d 名分\n', target_year, size(data, 5));
end

fprintf('\n=== 全年度の個人波形作成完了 ===\n');

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":25.2}
%---
