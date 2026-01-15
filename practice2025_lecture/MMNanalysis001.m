%[text] 入力
%[text] ord, path, prefix(年代指定)
%[text] やること
%[text] dataをhzごとに振り分ける
%[text] Triggerを抽出する
%[text] ordmemo
%[text] 16
%[text] 321
%[text] 231
%[text] 132
%[text] 
%[text] 17
%[text] 123
%[text] 213
%[text] 312
%[text] 231
%[text] 
%[text] 18
%[text] 312
%[text] 132
%[text] 
%[text] 22
%[text] 123
%[text] 213
%[text] 312
%[text] 
%[text] 23
%[text] 231
%[text] 321
%[text] 
%[text] 24
%[text] 321
%[text] 123
%[text] 231
%[text] 213
%%
% function [data0, trg] = processdata0(folderpath, ord, prefix)
% disp(folderpath)
% % データが置いてあるpath
% folder = folderpath
% % current directoryのファイル名を取得
% files = dir(fullfile(folder,'*.mat'));
% disp(files)
% % 複数の .mat ファイル名を取得
% names = {files.name}; % 名前の文字列をcell配列型にまとめる
% %lamda関数  sと言う引数をもち、1~8文字目まで出力。スライスを出力。
% % names = {'20xxxxxx_A', '19yyyyyy_B', '20zzzzzz_C', ...};
% 
% targetYear = prefix;                      % 抽出したい最初の2桁（3文字列）
% prefixes = cellfun(@(s) s(1:8), names, 'UniformOutput', false);
% 
% 
% fn0_filtered = unique(prefixes);
% nPrt = numel(fn0_filtered);
% 
% if nPrt == 0
%     error('prefix="%s" に一致する参加者IDが見つかりません。ファイル名先頭を確認してください。', targetYear);
% end
% 
% 
% 
% data0=zeros(31, 800, 1000, 3, length(fn0_filtered)); % (ch,time,epoch,freq,prt)の順番
% % ただし，prtは年度によって変わるので任意の数値が定義できるようにしておく
% data0(:,:,:,:,:) = 0;
% 
% % data0は(ch,time,epoch,freq,prt)の構造
% % EEG.dataは (chl,time,epoch)の構造
% % ord(part,ord_freq)
% for prt=1:size(fn0_filtered,2)
%     % eeglabで作成したデータを読み込む。
%     %[ ]は文字列連結できる印.
%     fn0_=dir(fullfile(folder,[char(fn0_filtered{prt}), '*.mat'])); % 各参加者の課題１～３を選択します
%     for dev =1:3
%         load(fullfile(folder, fn0_(dev).name));
%         data0(:, :, :, dev, prt)=EEG.data(:,:,1:1000); 
% 
%         % data0(＿＿,＿＿,＿＿,＿＿,＿＿)=EEG.data; %
%         % 本来なら右辺はEEG.dataで良いが2024年度3人目1回目の課題で上書きして1,000以上epochがあるため，EEG.data(:,:,1:1000);で1:1000で指定した
%     end
% end
% % まず入れ物を準備する
% trg=zeros(1000 ,3 ,size(fn0_filtered,2)); % (epoch,dev,prt)
% % epoch: 1 -> standard, 2-> deviant
% % dev: 1-> 1010 Hz, 2-> 1030 Hz,  3-> 1050 Hz
% % prt: 参加者
% 
% for prt = 1:size(fn0_filtered, 2) % participants
%     fn0_ = dir(fullfile(folder,[fn0_filtered{prt} '*.mat'])); % 課題】何をしているか日本語で説明せよ。
%     %Ans. データに対応する参加者IDを取り出し、それに対応するmatファイルのパスのリストを作成している
%     for dev = 1:3 % 課題順
%         load(fn0_(dev).name); % 課題毎にデータ読み込む
% 
%         % イベントタイプのセル配列を取得
%         event_types = {EEG.urevent.type};      
% 
%         % 'S  1'と'S  2'のインデックスを取得
%         % 標準刺激であるトリガーがどこにあるかのIndexたち(Indices)を取得している。
%         S1_indices = strcmp(event_types, 'S  1');
%         S2_indices = strcmp(event_types, 'S  2');
% 
%         % 各イベントに対するトリガー値を設定
%         %標準刺激か逸脱刺激かを示すトリガー値の箱を作成している。
%         trg_values = zeros(length(event_types), 1); % 【課題】何をしているか日本語で説明せよ。
% 
%         trg_values(S1_indices) = 1; % 【課題】何をしているか日本語で説明せよ。
%         %Ans 標準刺激が含まれるindexに標準刺激を示す1を代入している
%         trg_values(S2_indices) = 2; % 【課題】何をしているか日本語で説明せよ。
%         %Ans 逸脱刺激が含まれるindexに標準刺激を示す1を代入している
% 
%         % トリガー行列の更新
%         valid_triggers = trg_values(trg_values > 0);  % 【課題】何をしているか日本語で説明せよ。
%         %標準刺激、逸脱刺激の種類に関わらず、トリガーがある部分の行列を生成する。
%         trg(:, ord(prt, dev), prt) = valid_triggers(1:1000); 
%     end
% end
% save('2024data_mmn.mat','data0','trg',"-v7.3")
% end
% folderpath = '/Users/nyanziba/school/脳計測科学/MATLAB_EEG/Dataset/2024Autumn_mmn_set'
% ord = [3 2 1
%     1 2 3
%     2 3 1
%     2 1 3
%     ]
% 
% data02024=zeros(31, 800, 1000, 3, 4);
% data02024, trg = processdata0(folderpath, ord, "24")
%%
%[text] # 
%[text] # 個別データを整理し，一つのデータにまとめる。
%[text] EEGLABで作成したデータ（EEG.data）は，1条件ずつのデータで，
%[text] (chl,time,epoch)の3次元構造になっています。
%[text] サイズは(31, 800, 1000)です。
%[text] これを data0-\>(chl,time,epoch,freq,prt)に５次元構造に統合し，解析していきます。
%[text] ただし，
%[text]  freq（deviant条件）: 1010 Hz　-\> 1, 1030 Hz -\> 2, 1050 Hz -\> 3 ,
%[text] 例えば，計測順が25112801，25120501，25121201，25121202のとき，
%[text] prt（participants）: 25112801 -\> 1, 25120501 -\> 2, 25121201-\> 3, 25121202-\> 4とします。
%[text] ## 参加者の順番を定義する。
%[text] fn0に参加者のIDを定義する。
%[text] % 例：「2025年04月15日 1人目」は'25041501'を意味しました。
%[text] % この先頭 8文字の参加者IDをfn0に定義しましょう。
%[text] % fn0 = {'25', '04', '15', '01'};
%[text]                         ↑            ↑                 ↑                 ↑
%[text] % 2025年度の各参加者ID（頭8文字）を手入力で入れる
%[text] 
%[text] でもよいですが，
%[text] 少しずつ関数の使い方を学んでいきましょう。
%[text] **やりたいこと**
%[text]  **「参加者のIDはファイル名になっているので，複数の .mat ファイル名から、参加者番号（先頭8文字）だけを取り出し、重複を除いて一覧にする」**
% データが置いてあるpath
folder = '/Users/nyanziba/school/脳計測科学/MATLAB_EEG/Dataset/2016Autumn_mmn_set' %[output:89b60afa]
% current directoryのファイル名を取得
files = dir(fullfile(folder,'*.mat'));
% 複数の .mat ファイル名を取得
names = {files.name}; % 名前の文字列をcell配列型にまとめる
%lamda関数  sと言う引数をもち、1~8文字目まで出力。スライスを出力。
prefixes = cellfun(@(s) s(1:8), names, 'UniformOutput', false);
fn0 = unique(prefixes);
%%
%[text] ## data0の定義
%[text] 個別データで解析するよりもまとまっていた方が解析しやすいので，
%[text] データをまとめていきます。
%[text] EEG.dataは，1条件ずつのデータで，(chl,time,epoch)の3次元構造になっています。
%[text] サイズは(31, 800, 1000)です。
%[text] これを data0-\>(chl,time,epoch,freq,prt)に５次元構造に統合しましょう。
%[text] まず入れ物を準備します。
%[text] zerosでメモリの確保。
%[text] ＿＿に具体的な数字を入れてください。
data0=zeros(31, 800, 1000, 3, length(fn0)); % (ch,time,epoch,freq,prt)の順番
% ただし，prtは年度によって変わるので任意の数値が定義できるようにしておく
data0(:,:,:,:,:) = 0;
%%
%[text] ## 課題の順番を並べ替える
%[text] freq: 1 \<-\> 1010 Hz, 2 \<-\> 1030 Hz, 3 \<-\> 1050 Hz
%[text] 次に参加者ごとに課題を行った順番を定義します。
%[text] deviantが1010Hzの時を1010条件(1)，1030Hzの時を1030条件(2)，1050Hzの時を1050条件(3)とします。
% ord(part, ord_freq) の意味：
%   行：参加者
%   列：課題順
% （1 -> 1010 Hz, 2 -> 1030 Hz, 3 -> 1050 Hz条件とします ）
%
% 例：
%   1人目は 1010 → 1030 → 1050 Hz の順で課題を受けたなら，
%   ord(1,:) = [1 2 3] となります。
% 
% % -------2025年度
% 
% ord = [ 1 2 3
%         2 3 1
%         3 1 2
%         1 3 2
%         2 1 3
%         3 2 1
%         1 2 3]
% % --------2024年度
% ord = [
%     3 2 1
%     1 2 3
%     2 3 1
%     2 1 3
%     ]
% % -------2023年度
% ord = [2 3 1
%         3 2 1]
% % -------2022年度
% 参加者ごとの課題順を定義する
% ord = [1 2 3; 2 1 3; 3 1 2;]; % 2022年度の例
%2018年度
% ord = [3 1 2; 1 3 2;]
%2017年度
% ord = [1 2 3; 2 1 3; 3 1 2; 2 3 1]
%2016年度
ord = [3 2 1; 2 3 1; 1 3 2] %[output:5d2d9cd2]
% 自分たちの実験ノートを見て、
% 参加者ごとの課題順をここに正しく埋めること。
%%
%[text] ## dataを一つにまとめる
%[text] 参加者ごと，課題ごとにEEG.dataを読み込んで，
%[text] data0に代入していきます。
% data0は(ch,time,epoch,freq,prt)の構造
% EEG.dataは (chl,time,epoch)の構造
% ord(part,ord_freq)
for prt=1:size(fn0,2)
    % eeglabで作成したデータを読み込む。
    %[ ]は文字列連結できる印.
    fn0_=dir(fullfile(folder,[char(fn0(prt)), '*.mat'])); % 各参加者の課題１～３を選択します
    for dev =1:3
        load(fullfile(folder,fn0_(dev).name))
        data0(:, :, :, ord(prt, dev),prt)=EEG.data(:,:,1:1000); 
        
        % data0(＿＿,＿＿,＿＿,＿＿,＿＿)=EEG.data; %
        % 本来なら右辺はEEG.dataで良いが2024年度3人目1回目の課題で上書きして1,000以上epochがあるため，EEG.data(:,:,1:1000);で1:1000で指定した
    end
end
%[text] **for文を書いたときは**
%[text] 　prt=1 をまずコマンドウィンドウで実行し，1行ずつ走らせて問題なければループを回すこと。
%[text] **エラーが出たの対処法**
%[text] 小刻みに動かす：１行書いたら実行、を癖にする（最後に一括実行は事故る）
%[text] for分の場合は，for文の変数この場合はprtがどこまで進んだか確認，
%[text] どの様なエラーが出ているか確認し，関数の使い方，データサイズを確認すること。
%[text] 
%%
%[text] ## 関数の説明 strcmp -\> 文字列の比較
%[text] [`tf = strcmp(s1,s2)`](https://jp.mathworks.com/help/matlab/ref/strcmp.html#d126e1404145)
%[text] [`tf`](https://jp.mathworks.com/help/matlab/ref/strcmp.html#btwfvmr-tf) `= strcmp(`[`s1,s2`](https://jp.mathworks.com/help/matlab/ref/strcmp.html#btwfvmr-s1s2)`)` は、`s1` と `s2` を比較し、両者が同一の場合は `1` (`true`) を返し、そうでない場合は `0` (`false`) を返します。
s1 = 'Yes';
s2 = 'Yes';
tf = strcmp(s1,s2) %[output:4cc70105]

clear s1 s2 tf
%%
%[text] ## triggerの抽出（trgの定義）
%[text] 次に各課題の1000epochsで何の刺激が提示されたか分かっていないので調べる。
%[text] triggerの入っている場所 ，EEG.urevent.type を確認すると，
%[text] 2025年度はEEG.ureventは1001\*1のサイズで，boundary，S1，S2が入っています。
%[text] 2024年度以前はEEG.ureventは2000個近くある。
%[text] S  1 -\> standard，S  2 -\> deviant ，R1，R2はresponseの１と2のボタンです。
%[text] 解析の流れはppt参照すること。
%[text] 注意：「S  1」,と「S  2」の「S」と「数字」の間には半角スペースが2つ入っています。
% まず入れ物を準備する
trg=zeros(1000 ,3 ,size(fn0,2)); % (epoch,dev,prt)
% epoch: 1 -> standard, 2-> deviant
% dev: 1-> 1010 Hz, 2-> 1030 Hz,  3-> 1050 Hz
% prt: 参加者

for prt = 1:size(fn0, 2) % participants
    fn0_ = dir(fullfile(folder,[fn0{prt} '*.mat'])); % 課題】何をしているか日本語で説明せよ。
    %Ans. データに対応する参加者IDを取り出し、それに対応するmatファイルのパスのリストを作成している
    for dev = 1:3 % 課題順
        load(fullfile(folder,fn0_(dev).name)); % 課題毎にデータ読み込む
        
        % イベントタイプのセル配列を取得
        event_types = {EEG.urevent.type};      
        
        % 'S  1'と'S  2'のインデックスを取得
        % 標準刺激であるトリガーがどこにあるかのIndexたち(Indices)を取得している。
        S1_indices = strcmp(event_types, 'S  1');
        S2_indices = strcmp(event_types, 'S  2');
        
        % 各イベントに対するトリガー値を設定
        %標準刺激か逸脱刺激かを示すトリガー値の箱を作成している。
        trg_values = zeros(length(event_types), 1); % 【課題】何をしているか日本語で説明せよ。
       
        trg_values(S1_indices) = 1; % 【課題】何をしているか日本語で説明せよ。
        %Ans 標準刺激が含まれるindexに標準刺激を示す1を代入している
        trg_values(S2_indices) = 2; % 【課題】何をしているか日本語で説明せよ。
        %Ans 逸脱刺激が含まれるindexに標準刺激を示す1を代入している
        
        % トリガー行列の更新
        valid_triggers = trg_values(trg_values > 0);  % 【課題】何をしているか日本語で説明せよ。
        %標準刺激、逸脱刺激の種類に関わらず、トリガーがある部分の行列を生成する。
        trg(:, ord(prt, dev), prt) = valid_triggers(1:1000); 
    end
end
clear prt dev event_types S1_indices S2_indices trg_values valid_triggers
%[text] 
%%
%[text] ## workspaceの保存
%[text] save　ワークスペース変数をファイルに保存 save(filename,variables)
% 出力先フォルダ設定（ユーザーが `outputDir_mmn` を事前に定義していなければ既定値を使用）
% 既定値: プロジェクト直下の "MMN_raw" フォルダ
if ~exist('outputDir_mmn','var') || isempty(outputDir_mmn)
    thisFileDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisFileDir)); % .../脳計測科学
    outputDir_mmn = fullfile(projectRoot, 'MMN_raw');
end
if ~exist(outputDir_mmn, 'dir')
    mkdir(outputDir_mmn);
end

% 年度に応じた出力ファイル名を自動決定（folder パスから 20xx を抽出）
yrToken = regexp(folder, '(20\d{2})', 'match', 'once');
if isempty(yrToken)
    % フォルダ名から取得できない場合は固定名で保存
    outFile = 'data_mmn.mat';
else
    outFile = sprintf('%sdata_mmn.mat', yrToken);
end

% フルパスで保存
save(fullfile(outputDir_mmn, outFile), 'data0', 'trg', '-v7.3')

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":29.6}
%---
%[output:89b60afa]
%   data: {"dataType":"textualVariable","outputData":{"name":"folder","value":"'\/Users\/nyanziba\/school\/脳計測科学\/MATLAB_EEG\/Dataset\/2016Autumn_mmn_set'"}}
%---
%[output:5d2d9cd2]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ord","rows":3,"type":"double","value":[["3","2","1"],["2","3","1"],["1","3","2"]]}}
%---
%[output:4cc70105]
%   data: {"dataType":"textualVariable","outputData":{"header":"logical","name":"tf","value":"   1\n"}}
%---
