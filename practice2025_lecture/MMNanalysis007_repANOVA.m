%[text] # 統計解析(t検定～Ⅰ要因３水準の反復測定ANOVA)
%[text] 目標：MMNが起きているチャンネルを探し，条件間で比較する。
%[text] １．データを扱いやすいように変形させる。
%[text] data\_aveは，３次元目にstd条件とdev条件があるので，
%[text] (dev条件) - (std条件)とデータを変形する。
%[text] ２．さらに，このときの波形を観察してどの期間にMMNが発生してそうか，「MMNが出ていそうな時間窓」を決める。
%[text] ３．その時間窓の平均振幅が 0 と異なるか評価する（1標本t検定）。
%[text] ４．有意だったチャンネルで deviant 3条件（1010/1030/1050 Hz）を比較（反復測定ANOVA）。
%%
%[text] ## 0. データの準備
% TODO: EEG電極座標とチャンネル名を読み込む
load('EEG_locs_new.mat');     % sloc, ch_name が入っている想定

% TODO: data_ave を workspace に用意（既にあるならOK）
% whos data_ave

%%
%[text] ## １．データを扱いやすいように変形させる。
%[text] MMNの大きさを条件毎に比較したいのでdeviant条件からstandard条件を引く。
%[text] 変数は，data\_av % (chl, time, \[std or var\],\[1010 or 1030 , 1050\],prt)
%[text] を使えばよい。
% TODO: std/dev の次元（3次元目）を確認して差分を取る
% data_dif = squeeze( ?????? );  % (ch, time, freq, participant) を目標


% チェック
disp('size(data_dif) ='); disp(size(data_dif));
% 期待： [nCh, N, 3, nP]
%%
%[text] ## ２．さらに，このときの波形を観察してどの期間にMMNが発生してそうか，時間帯を調べる。
%[text] チャンネル毎にMMN（deviant - standard）の平均波形を書いてください。006.mlxで書いた平均波形と一致することを確認してください。
fs = 1000; % サンプリング周波数

N=size(data_dif,2); % データ長 data(ch,time,blk,knd);% 1->1050 Hz, 2->1010 Hz

marker=zeros(2,2);
marker(:,1)=0.0;
marker(1,2)=-20;
marker(2,2)=20;
 
sz = get(0, 'ScreenSize');

x=linspace(-0.1,(N-100)/fs,N);

figure('Position', [sz(1) sz(2) sz(3) sz(4)])
% figure
for chl=1:size(sloc,1)
    subplot('position',[sloc(chl,1) sloc(chl,2) 0.07 0.07])
    hold on

    % TODO: 平均の次元（participant）を確認しながら mean を書く
    % y1010 = mean( ?????? , ??? , 'omitnan');
    % y1030 = mean( ?????? , ??? , 'omitnan');
    % y1050 = mean( ?????? , ??? , 'omitnan');

    plot(x, y1010, 'Color', [255  75   0]/255, 'LineWidth',2);
    plot(x, y1030, 'Color', [  0  90 255]/255, 'LineWidth',2);
    plot(x, y1050, 'Color', [  3 175 122]/255, 'LineWidth',2);
    axis('ij')

    ymx=10;
    title([ ch_name{chl} '(' num2str(chl) ')' ]) % titleを付けます。
%     title(ch)
    xlim([-0.1 0.7])
    ylim([-ymx ymx])
    plot(marker(:,1),marker(:,2),'k:','LineWidth',2);
    set(gca,'box','on');
    set(gca,'linewidth',1)
end
legend( '1010Hz',  '1030Hz', '1050Hz')
makeSubplotsClickable();
%%
%[text] 全チャンネル表示からめぼしいチャンネルを選択してfigureを描いてください。
% TODO: 候補チャンネル番号を自分で入れる（複数OK）
chl_list = [5,14];  % 例

figure('Position', [sz(1) sz(2) sz(3) sz(4)]) %[output:479e4244]
hold on %[output:479e4244]
for chl=chl_list
    % TODO: ここも participant 平均の取り方を再確認
    plot(x, mean(data_dif(chl,:,1,:),4,'omitnan'), 'Color', [255  75   0]/255, 'LineWidth',2); %[output:479e4244]
    plot(x, mean(data_dif(chl,:,2,:),4,'omitnan'), 'Color', [  0  90 255]/255, 'LineWidth',2);
    plot(x, mean(data_dif(chl,:,3,:),4,'omitnan'), 'Color', [  3 175 122]/255, 'LineWidth',2);
    axis('ij')
    title([ch_name{chl} '(' num2str(chl) ')'])
end
    ymx=7;
    title([ ch_name{chl} '(' num2str(chl) ')' ]) % titleを付けます。 %[output:479e4244]
%     title(ch)
    xlim([-0.1 0.7]) %[output:479e4244]
    ylim([-ymx ymx]) %[output:479e4244]
    plot(marker(:,1),marker(:,2),'k:','LineWidth',2); %[output:479e4244]
    set(gca,'box','on'); %[output:479e4244]
    set(gca,'linewidth',1) %[output:479e4244]

legend( '1010Hz',  '1030Hz', '1050Hz') %[output:479e4244]
%%
%[text] % TODO: 時間窓を自分で決める（例：0.150〜0.200 s）
%[text] 秒→data pointに変更する。(+100points)
% 目標：mmn_amp (ch, freq, participant)

% size mmn_amp 今年は(31,3,20)です， (ch,freq,part)
%%
%[text] ## ３．その期間のamplitudeが０ではないことを統計学的に評価する（t検定）。
%[text] 
%[text] #### `関数    ttest`
%[text] `1 標本およびペア標本 t 検定`
%[text] [`[h,p,ci,stats] =` ](https://jp.mathworks.com/help/stats/ttest.html#d122e704342)[<u>`ttest(x)`</u>](https://jp.mathworks.com/help/stats/ttest.html#d122e704193)
%[text] `平均がゼロで分散が未知の正規分布であるという帰無仮説の検定の判定を返します。`
%[text] `p：p値`
%[text] `規定値（デフォルト）の有意水準は，0.05`
%[text] `有意の時 h は 1、それ以外は 0`
%[text] `信頼区間 ci` 
%[text] `statsの中身，`
%[text] `tstat — 検定統計量の値。（stats.tstatで取り出せる，以下同様）`
%[text] `df — 検定に対する自由度。`
%[text] `sd — 母標準偏差の推定値。`
%[text] 
%[text] ただし，ttest は「1次元目」に対して検定するので、 (ch,freq,participant) → (participant,ch,freq) に変換が必要。
%[text] 
%[text] #### `関数    permute`
%[text] `permuteを用いて，mmn_ampの次元の並べ替えをします。`
%[text] `関数ttestは1次元目を評価してくれるので，partを３ー＞１次元目に移動`
%[text] `(chl, freq, part) => (part, chl, freq)`
%[text] 
%[text] [`B = permute(A,dimorder)`](https://jp.mathworks.com/help/matlab/ref/permute.html#f93-999230) `　` 
%[text] `今回は、ttestで比較したい３次元目を1次元目に持ってきたいので，`
% TODO: (ch,freq,participant) → (participant,ch,freq) に並べ替える
A = permute(mmn_amp, [_____]);
%[text] 
[~,p,~,stats] = ttest(A); 
%%
%[text] ### t検定の結果を表示（t-map）
%[text] 1010, 1030 , 1050条件ごとに，チャンネル位置でｔ値をカラーマップ（グラデーション）で表示します，
%[text] ただし，ｐ＜0.05の有意のチャンネルだけ表示，ｐ＞0.05のチャンネルはｔ値を０とします）

amp=10;
thr=0.05;

figure; hold on
for f = 1:3
    subplot(3,1,f); hold on

    % TODO: t値マップを取り出す（stats.tstat の次元を見て書く）
    % Map = ??????

    % TODO: p>thr を 0 にする（同じ形状になるように）
    % Map( ?????? ) = 0;

    drawPatch34rwb(sloc, Map(:,:), 0.03, [-amp amp]);
    colorbar('Ticks',[0,0.5,1],'TickLabels',{-amp 0 amp})
    axis tight; axis off;

    if f==1, title('tmap(1010Hz)'); end
    if f==2, title('tmap(1030Hz)'); end
    if f==3, title('tmap(1050Hz)'); end
end
%[text] 
%%
%[text] ## ４．t検定で有意になったチャンネルにおいて，３条件（deviant）をanovaで比較する。
%[text] ## 4.1. データの様子（分布）を見てみる。
% TODO: 有意チャンネルを自分で決める（例：t-mapから）
chl=14 %[control:slider:9c63]{"position":[5,7]}

% 条件（freq）×参加者の行列（nP × 3）を作る
X = ＿＿＿＿＿＿;   % (participant, freq)

% ---- 4.1 分布確認（散布図）
nP = size(X,1);

figure
hold on
scatter(1:nP, X(:,1), [],'MarkerEdgeColor', [255  75   0]/255,'Marker','o')
scatter(1:nP, X(:,2), [],'MarkerEdgeColor', [  0  90 255]/255,'Marker','^')
scatter(1:nP, X(:,3), [],'MarkerEdgeColor', [  3 175 122]/255,'Marker','s')
axis('ij')
xlabel('participants (number)')
ylabel('amplitude (microV)')
xlim([0, 22])
legend('1010','1030','1050','Location','best')
%%
%[text] ## ４.2. バーグラフを書く（平均・標準誤差）。
%[text] ## 条件毎の電圧の平均値のバーグラフ
%[text] `エラーバーを標準誤差にします。`
%[text] `関数 std    標準偏差　バラツキ具合を表す。`
%[text] [`S = std(A,w,dim)`](https://jp.mathworks.com/help/matlab/ref/std.html?s_tid=srchtitle#d122e1163169)
%[text] `w：重み付けスキームを指定します。w = 0 (既定値) の場合、S は N-1 で正規化されます。w = 1 の場合、S は観測値の数 N で正規化されます。`
%[text] `dim：次元に沿った標準偏差を返します。`
%[text] `標準誤差（standard error）= std / √n`
%[text] 
% mmn_amp (31,3,○○) (ch,freq,part)
chl=14 %[control:slider:768c]{"position":[5,7]}

sem = std(mmn_amp(chl,:,:),0,3,'omitnan') ./ sqrt(sum(~isnan(mmn_amp(chl,:,:)),3));

figure
bar(1:3,mean(mmn_amp(chl,:,:),3,'omitnan'))
hold on
er = errorbar(1:3,mean(mmn_amp(chl,:,:),3,'omitnan'),sem,sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
set(gca,'XTick',[1 2 3]); % xのそれぞれの値の所に
set(gca,'XTickLabel',{'1010', '1030', '1050'}); % 実際の値を定義する。
xlabel('deviant condition (Hz)')
ylabel('amplitude (microV)')
axis('ij')
%%
%[text] ## ４.2. violinplotを書く。
%[text] 散布図と最大値，最小値と「25パーセンタイル（第一四分位数）」、「50パーセンタイル（中央値）」、「75パーセンタイル（第三四分位数）」を表示する。
%[text] 
chl=14 %[control:slider:1c7f]{"position":[5,7]}
figure
y = squeeze(mmn_amp(chl,:,:))';
vs = violinplot(y);
hold on
s = swarmchart(1:1:3,y,'filled');
% set(s, 'SizeData', 15); % 例：20 くらい?
set(s, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
boxchart(y)
set(gca,'XTickLabel',{'1010', '1030', '1050'}); % 実際の値を定義する。
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('deviant condition (Hz)')
ylabel('amplitude (\muV)')
axis('ij')
% xlim([0.5, 3.5]);
%%
%[text] ## 4.4. 3条件の統計学的に比較する。
%[text] 反復測定の分散分析を使用します。
% chlの定義
chl=5; % 適当

% 条件（周波数）×参加者のデータを作成
X=squeeze(mmn_amp(chl,:,:))';

% 参加者間情報を格納する変数を定義
% fn0(paticipantsの名前)を作成
% fn0を既に作成していればそのまま使用する。
% もし，作成していなければ適当な名前で良いので，以下をアンコメント化して使う。
for prt=1:size(X,1)
    fn0{prt}={['prt' num2str(prt)]};
end
% 参加者間要因のデータをテーブル配列形式で作成
between = table(fn0',X(:,1),X(:,2),X(:,3),...
    'VariableNames',{'prt','Hz1010','Hz1030','Hz1050'});

% 参加者内要因をテーブル配列形式で作成
within = table({'Hz1010';'Hz1030';'Hz1050'},'VariableNames',{'Deviant'});

% 反復予測モデルを近似
rm = fitrm(between,'Hz1010-Hz1050~1','WithinDesign',within);

%[text] ### 球面性の検定
%[text] Mauchlyの球面性検定
mau = mauchly(rm)
%[text] #### **解釈**
%[text] - `mau.pValue > .05` → 球面性OK（通常の自由度）。
%[text] - `mau.pValue < .05` → 球面性NG（GGやHFで自由度補正が自動で表示される列を参照）。
%[text] - Greenhouse-Geisser（クリーンハウス・ゲイザー）やHuynh-Feldt（ホイン・フェルト）による自由度の調整が必要となり，調整した値を用いて*p*値を求める必要になる。
%[text] - 基本的にはGreenhouse-Geisser（クリーンハウス・ゲイザー）で良いです。 \
%[text] #### 結果の書き方の例
%[text] 反復測定分散分析の結果、条件間の主効果がみられた
%[text] （*F*(2, 64) = 20.43, *p* \< .001）。
%[text] 
%[text] #### もし Mauchly 検定が有意（p \< .05）だった場合
%[text] 反復測定 ANOVA のイプシロン調整
%[text] イプシロン係数（Greenhouse-Geisser）を自由度と誤差それぞれにかけて，報告する必要があります。
ep = epsilon(rm)
%[text] 球面性の仮定が棄却されたため、Greenhouse–Geisser 補正を適用した。
%[text] Greenhouse–Geisser の ε = .88 で補正した結果、条件間の主効果があった
%[text] （*F*(1.76, 56.32) = 20.43, *p* \< .001）。
%[text] 
%[text] ### 反復測定の分散分析を実行
ranovatbl = ranova(rm)

%%
%[text] # 効果量
%[text] 
%[text] ANOVA系でよく使われる効果量（effect size）は主に以下の2種類です。
%[text:table]{"columnWidths":[-1,-1,139,-1]}
%[text] | **名称** | **記号** | **定義** | **一般的な解釈** |
%[text] | --- | --- | --- | --- |
%[text] | **η² (eta squared)** | η² = SS\_effect / SS\_total | 全体の分散に対して要因が説明する割合 | R²に似た概念 |
%[text] | **偏η² (partial eta squared)** | ηp² = SS\_effect / (SS\_effect + SS\_error) | 要因が誤差を除いた分散のうちどの程度を説明しているか | 🔸**通常はこちらを報告** |
%[text:table]
%[text] 
%[text] 反復測定ANOVAでは，被験者内要因の「誤差（Error）」が明示的にあるため，
%[text] **偏η² (partial eta squared)** の報告が標準的です。
%[text] 
%[text] ## MATLABの出力から計算する
%[text] `ranova(rm)` の結果には `SumSq`（平方和）が含まれています。
%[text] これを使って偏η²を計算します。

SS_effect = ranovatbl.SumSq(1);  % 要因（frequency）
SS_error  = ranovatbl.SumSq(2);  % 誤差項


eta_p2 = SS_effect / (SS_effect + SS_error);
fprintf('Partial eta-squared = %.3f\n', eta_p2);
%[text] 
%[text] ## 効果量の目安（Cohen, 1988）
%[text:table]
%[text] | **効果量** | **ηp²の範囲** | **解釈** |
%[text] | --- | --- | --- |
%[text] | 小 (small) | 0\.01〜0.06 | 弱い効果 |
%[text] | 中 (medium) | 0\.06〜0.14 | 中程度 |
%[text] | 大 (large) | \> 0.14 | 強い効果 |
%[text:table]
%[text] 
%[text] ## レポートでの書き方例
%[text] 今回の解析結果とは異なりますので，そのまま記述しない様に気を付けてください。
%[text] 
%[text] 反復測定分散分析の結果、条件の主効果は有意であった
%[text] （*F*(1.76, 56.32) = 20.43, *p* \< .001, ηp² = .39, Greenhouse–Geisser 補正）。
%[text] （注意）ηp²のpは下付き，²はηの上付きです。
%[text] 
%[text] 反復測定分散分析の結果、**条件の主効果は有意ではなかった**
%[text] （*F*(1.76, 56.32) = 1.42, *p* = .25, ηp² = .04, Greenhouse–Geisser 補正）。
%[text] 
%%
%[text] ## 事後検定（多重比較）
%[text] 
% 'ComparisonType' ? 使用する棄却限界値の種類
% 'tukey-kramer' (既定値) | 'dunn-sidak' | 'bonferroni' | 'scheffe' | 'lsd'
T0 = multcompare(rm,'Deviant')
T1 = multcompare(rm,'Deviant','ComparisonType','dunn-sidak')
T2 = multcompare(rm,'Deviant','ComparisonType','bonferroni')
T3 = multcompare(rm,'Deviant','ComparisonType','Scheffe')
T4 = multcompare(rm,'Deviant','ComparisonType','lsd')

%%
%[text] ## バイオリンプロットとsigstar

figure

y = squeeze(mmn_amp(chl,:,:))';

vs = violinplot(y);
hold on
s = swarmchart(1:1:3,y,'filled');
% set(s, 'SizeData', 15); % 例：20 くらい?
set(s, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
boxchart(y)

axis('ij')   % ← ここで下向き表示に固定

% sigstar
groups = {[1 2],[1 3],[2 3]};
pvals  = [T0.pValue(1), T0.pValue(2), T0.pValue(4)]; % pvalsは'ComparisonType'を自分で選ぶ

ss = sigstar(groups, pvals, 0, 'top');   % 括弧と星
% ss = sigstar(groups, pvals, 0, 'auto');  % 自動（画面上）
% ss = sigstar(groups, pvals, 0, 'bottom');

set(ss(:,1), 'Color', 'k', 'LineWidth', 1);
set(ss(:,2), 'FontSize', 12);      % 星の大きさ

% Y=get(ss(1,1),'YData'); Y(1)=-7.5; set(ss(1,1),'YData',Y); % 有意差が出た場合の足の長さを調整します。


% 軸・ラベル
set(gca,'XTickLabel',{'1010', '1030', '1050'});
% set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('deviant condition (Hz)')
ylabel('amplitude (\muV)')









%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":29.6}
%---
%[control:slider:9c63]
%   data: {"defaultValue":14,"label":"chl","max":31,"min":1,"run":"Section","runOn":"ValueChanging","step":1}
%---
%[control:slider:768c]
%   data: {"defaultValue":14,"label":"chl","max":31,"min":1,"run":"Section","runOn":"ValueChanging","step":1}
%---
%[control:slider:1c7f]
%   data: {"defaultValue":14,"label":"chl","max":31,"min":1,"run":"Section","runOn":"ValueChanging","step":1}
%---
%[output:479e4244]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAecAAAESCAYAAADHUh3+AAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl8FdX1\/hKorWCrWBEVJAsIKJuiFFEkCQEUxAWpGy4sAZdWBStaRZBFUVvccF8IEBeUWhQ3ZAuERQWXIFp3TQJiXbBG\/wpq1eTfMy\/z3sy8OzN33pt5b+57Z34\/fkByl3O\/c+d+c+4995ycxsbGRvDDCDACjAAjwAgwAqFBIIfJOTS6YEEYgcAQKCwsRE1NjW\/t+92eb4JxQ4xAhiDA5JwhiuRhMAJ2CARFpEG1y5pkBBgBgMmZZwEjEFIEiPyMTyKWrxOBeiFXXRarDF7aCCnMLBYjEEoEmJxDqRYWKpsR8JMI7cjTrg8R7noboraYnLN5pvLYg0SAyTlIdLltRsAjAm5k5\/Z7Y3duxExlZaxxJ3KmNrzI5BEOLs4IZC0CTM5Zq3oeeBgRkCE6mTIypCnTjrGME9nLkHwY8WaZGIGwIsDkHFbNsFxZh4AMWcqQrg6cW3tef8\/knHVTkgecRgSYnNMIPnfNCMhsQ9uhZHUYM5YjS9Yr+Vr7sdZncub5ygikDgEm59RhzT0xAo4IuJGpLHxuZ8QylrUX5y+\/5JYdH5djBLIBASbnbNAyj1EJBPwgOZkzYjdy9moh+yG3EgpiIRmBFCLA5JxCsLkrRsAJAVmSkyVPt\/ZkvLnt5DU6gLn1w1pnBBgB7wgwOXvHjGswAoEhIEN0slvObm25\/d46SNmPgsDA4YYZgSxCgMk5i5TNQw0\/Am6EKUvMNNJE2krEsnfrJ\/yos4SMQPgQYHIOn05YoixHQPfClgmV6UaMTr93qytjOXttI8tVy8NnBKQRYHKWhooLMgKpRcB6VUpE1nYS6WWZnFOrM+6NEfALASZnv5DkdhiBkCIQlHUbVLshhZHFYgRSigCTc0rh5s4YgfQg4DeR+t1eelDhXhmB8CLA5Bxe3bBkjAAjwAgwAlmKQNaS8\/bt27NU5TxsRoARYASyA4F27dopO9CsJGci5ssvvxybNm1STnG1tbUmmQsKCpQbAwvMCDACjEAqEOjTpw9mz54NFUk6K8l548aNGDlypKa0tm3bpmKO+NZH3759TW299NJLvrWdiobog2jOnDlKYs+yp2KGxPfBuDPuiSCgz5t169YxOScCYDrq6OSsotJycnJMkDU2NqYDwoT7VBl7lj1htSdVkXFPCr6EK6uMOw1adfmz2nJmck74vU24Ih0pLF68GBMmTEi4jXRVZNnTgzzjzrgnggCTcyKopbmO6kpLM3zcPSPACDACoUdA9XWeLWeFvflC\/3awgIwAI8AIpAkBJuc0AZ9Mt6orLZmxc11GgBFgBLIBAdXXebacFbOcVXcIy4ZFgcfICDAC6UeAyTn9OvAsgcpKY3L2rG6uwAgwAlmIgMrrPKmLLWe2nLPwteUhMwKMQKYjwOSsoIZVVhpbzgpOOBY5UAQo2h9dz+MnsxHwevVV5XU+Iy1nYw5ca\/5bfeqqrDQm58xegHh03hFgcvaOmYo1mJxV1FqTzNY0dnZp7VQmZ4XVw6IzAoEgoJMzheMdMWJEIH1wo+lDoH\/\/\/qBANEzO6dNBUj2LiJjJOSlIuTIjoAQCTM5KqClhIZmcE4YuHBW9JH9nyzkcOmMpGAE\/EGBy9gPF8LbB5Bxe3UhJppOzlzPnhQsXRrOVqJhSTAoYLsQIZDgCTM6ZrWAv5Ezb3\/pD\/6bsg163w8OCpjJXqYykawWPHL\/03xudwNy2tY3tUCIGFZIxsENYWF4dliMsCDA5y2nCaXfRyahJpJ5bHTtnXdFIvJAzpaOlP8aHyVlufgRWKpEzZ2M+Z7KcVbCemZwDm0LcsKIIMDm7K05kvOi1nBxpk6lnR8BejiBJRi\/kTNaybj1zPmf3eZGSEomQs4pfVEzOKZlO3IlCCDA5OysrEavYulNpJVo7gjUeL6aDnI1IqO5bpMy2tttaweTshhD\/nhHITAQcyXnhdGBVRWYOXDSqgaOAkdOF4\/WyRhrLJlLPbVtbJKAdmXuxnJmcQzrV+Z5zSBXDYjECASLgSM63jgEqFwTYe8iaHjktNOTshIysJe51W5vJOWTz0SiOF29tFbe1Qww9i8YIpAUBJmcD7CEiZ9ltbbczaLac0\/JapadT1c8i0oMa98oIhBMBJmd1ydmNmNlyDuc7F5hUKpMzO4QFNi24YUURcCTnz+sUHVWCYrfJt62YyNkxNZZIPbczZ\/36q8yVKracE5wLKlZjclZRaywzIyBGgL215WaGm4e13oqs747bFSynbW3qS4aY2XKW023GlGJyzhhV8kAYATA5y00CN2tWb8WLs5adj49sX0bJ2VvbrMeMuUolNz0jpZicvaDFZRmBcCPA5Bxu\/SQrHW9rJ4ugQvWZnBVSFovKCLggwOSc2VOEyTmz9WsancrknEVq4qEyAlIIMDlLwaRsISZnZVXnXXAmZ++YcQ1GIKwIMDmHVTP+yMXk7A+OSrTC5KyEmlhIRkAKASZnKZiULcTkrKzqvAvO5OwdM67BCIQVASZnOc3IelBb0+4aWxd5cuu\/l0nXS2VlAo8Y+2RyltNvRpRSmZw5CElGTEEehI8IMDm7g5lI6ke34CPJ3HOWveNMI2NydtdvxpRgcs4YVfJAGAG+5+wyB5JJGWl359ktoImblc7k7P7i8j3ndu3cUQpRCbacQ6QMFiUUCDhZzmPuAxasC4WYKRFi2inA9BHirtwsYWOtRMlXNp+zSEIOQmJGhcmZyTkliwZ3wggEhQCTcwzZIMhZZHnLkLeTvr1EIeNt7aDenBC2y9vaIVQKi8QIJIgAk3Ow5CyyqGXImVNGJjihm6qx5ayY5Zycurk2I5B5CDA5q0vOMp7bbDln3jtrOyKVLecsUhMPlRGQQoC9taVg8pz6kVoN0iFMhphJBiZnOf1mRCmVyZkdwjJiCvIgfESAyVkOTLetaL0VJ6cuYxvJXKUSEb\/dKJic5fSbEaWYnDNCjTwIRkBDgMlZbiK4XW\/SW0lFEBKRxOytbUaFz5wVO3Nmy1luIeJS2YMAk3Nm65ot58zWr2l0bDlnkbJ5qBmPAJNzZquYyTkD9Wu3jcPknIHK5iFlLQJMzpmteibnDNOvUyxZlck5w9TEw2EEkkaAyTlpCEPdAJNzqNXjTTi3EHJMzt7w5NKMQJgRYHIOs3aSl43JOXkMQ9GCk6u\/LiCTcyhUxUIwAr4gwOQsB6NXb21j2E5jD3be3JwyUk4PsqUyzlvbCzkvXLgQ7Zq8tfW\/ZYHjcowAIxAOBJic3fXgV8pIY0\/J3HMOKivV9u3boyLSv0eOHIl169ZF13l3pMJTQhlytvuKIyh1RTtNFiPkuuVs\/NmECRNAf8L+8FWqsGuI5Us1Ak7kXLerPtXipLW\/\/Bat4vr3M2Wk3rhbQBM3Kz0ocp4zZw7oj\/Fhck7rlISn0HQ6Oc+ePRtt27bVJCfLWQXrmck5zRONuw8dAk7kPLZ6MR7cVh06mYMS6JouA3BNl1Jh815TRlobkdm2dvP3oTZltsuNfXs5cyZrWbeeN23apBE1k3NQs02yXRnLWm9K5TNnJmfJCcHFsgYBJueYqv0mZztClrGcnSYgp4x0fz2V2dZ2H0p8Cb7nnAhqXIcRUAsBJudgyFm09exmGbv9XrecZaxwfVReLGfjzFXZCKNxMDkrFr5TrWWTpWUEgkfAiZxnvluJim2bgxciJD2Man+Yb9vaqSBnmcxUTM4hmVx+ipGJlrOf+HBbjEAmIMDe2nJaTOTMmVNGymEbRKmMtpztAFN9uyOIiaBqm3U7gKp3gIp1QN2XsVF02BeYMhzI3wfIb63q6FhuGQSYnGVQEjvNOm0zu91+SeYqFfUr67HNlrOcfjOilMrkzA5hsSlIxFwyC6C\/nR4i5\/nnA8UHZ8T05UFYEGBylpsSbteb9FZkgozoZe2uacn2ZZScU0aa9ciWs2JnzkzOkQm8YB0w5j65RclY6uADgKVXsDXtHbnw1mByDq9u\/JCMLWc\/UFSkDbacFVGUjZhkKRdMTG4MZE2vuZpJOjkUw1GbyTkceghKCibnoJANYbtMziFUiqRIdlvZtG1NhEtnzPTQ+fP1TwEr33RueORRwCN\/luyci4USASbnUKrFN6GYnH2DMvwNMTmHX0d2EpZcF3EAMz5EzKP7i2voDmMznrA\/myZSr71NXUyyXXIm58yeAUzOma1f0+hUJucsUpNwqDlnyROztQEi6jnLgNuWiVFcM4WdxlScXzo5qyg7yyyPgNcwnKqv8+wQpphDmPxUzqySRKzkAGa0msn7mgjV6yNqS29j+ZXA4O5eW+Ty6USAyTmd6Keubybn1GGdtp5U\/6JKG3Bp7JhImba00bIeOa1rkZNfrf1NT16LVpjapQTF+xRClJXHTuzpiwHa7rY+lZOBAV3TOFjumhFgBJJGQPV1ni1ntpyTfgmCbiDqBLarHs2KyjWCtnvW9z8fffduLy1SlPQtNegM2rfgJZ\/XAW3ygZXzgPpPgT3bAAd0BPbNj\/ycH0aAEfAdASZn3yENvkGVlZaN95y1O80PuxOzPnPIkq7sVyZtRdtdzUr6DJpI+c0q4LYx9pP60IHAdSuDn\/Tcg3II6LmovewGKTfIAAVWeZ0nWNhyVsxyzjZy1qzmW+uxrWBxdBvbSMJbd9lb0ZX9xqFonwKp1z8QC3psAfBFnVT\/uGEN0L1YriyXyjgEiIi37voa1767GlVf1tiOjz48C1q0wrntD0O\/3+cjNycHVPfnhgbt381ycqXnfMaBaBkQk7OCGlZZadlGzmQ1j11ejdzei00zjUi3vNcIzTpe+2UtyqoXa4uU9aEyo9r3kpqlIoJOyOls1QJna9lOmj4nAVOXSMnKhcKNgE62lBWL5qfxIYI1Pk4fmMmMkt6NVR52kJLpK4x1VV7n2XJetw7t2HIO43sVlYms5vUFN5lkPHH\/g\/FEn7Pj5J639VWct\/nJuJ97saBFYUHnnQeMKZKAibaxq5cDd10QX5jOly+4A8jrFtnqXjhDbFVTuXnmxVyiZy4SIgQqtlVrH4theezel7DIF5QcTM5BIRtgu6orLUBoQtU0bWl3eLIcOfvGtvnczpPtFkYvBD32fmD+WjMUrg5iRMxXlYgJd\/pzwBFD47GlOmWCbXcm6FDNQy\/CdFhxE\/ywhGme57fYK9p13a6vk2qX2pvXa0RWbXmrvs7zmbNilrOXhUL1siX31WJ9m7mmYchsU9M2YukGcz1q5MPBk6SdxCh2tzXblS1B0za2yBImkr10vvNZMhE0OYyRNW18mKCVmr60jd1xhXmHx8sAiDxP2K8Lhh\/QFXkt9rKdp\/p2ObVtdRhb+2UNjmjVDhdvecaWyMfkHY4HDjvFi2jKlmVyVlB1qitNQcg9i6xZzUvmmpzA6JyZLGCZJ1mCFnlwm86fiVTpEREr\/ZzIlZy8ZK9K3XIusPohJmgZ5YasjN1cs7NWrb4RQXhjUx\/jNz+BNTvinctIrmu6DJD2xQgZ3NLiqL7Os+WsmOWcLQ5hf1lei9u\/N1u\/XramdctCZM3ItiM6fz6j21d49ICbgMdvsF8kDhsMXLtcehGJFrQjaC8k771XrpEEAnYWM1mnZKWm+7H7cCC5Mp2kmZzTPfsS6F9lpWUDOXu1mo3bz9bAIXaLkyxBn3gz8Ey1eZKV\/lyJVTsHRn9Yl5uPqubFqPjVKO1vekiOQd2A++UM\/VgHZInTNrnx8WqFJ\/BOcJXEEGi+5Oq4irJzK7EevdeiDwhyULN6jVtbyrRzaZXXedINW85sOXt\/2wOuMWN9La79j73VTGRMKSFnLI7PUKUTI2Wqom1oekTWDS1EHw2e5DwS2rre+BRKFvaMkq6xQvHPkXNinZDtGvMczGTGMOCV5+Kb47vQAc88+eZpTpVuKI872w0bMRtH5MWLfEibTnim7yh5QEJYksk5hEpxE0llpWW65Szy0DaeNTslrRDpnSzYTvsBk06rx3GVi5HT4utosaLWBagcPEJz\/IoL1WnwvibLeMzu811J2GneXXMKMGOE28w0\/P7ei4Fn72SC9gBZqoraWaIyzopWGbUPzR2RGO\/0wUlPQWtg9dWRn8uEkLUmg3HCgWSv2fkVBr8wTwqugft2xL2HniztSCnVaIoKqbzOE0QZZzkXFhZqqq+psY+yo7LSMp2cK16rR9nHZq9X3RqJxtjeEczbTQuhZnHvLb7iNKTlUixrPsRz56P\/uwDrmvdH1fhaHNjwMTBwtFwbz94F3HtRfFne5pbDL6BSF295GvfUbjK1nggx3\/A0MHmRs5A0J9c07ZzrRE3vQeVbwMMb7HeO7h0LHNvDuW0i6mc+fQeXvinYpbFUHdmuJ8rye5uuYlF9cjh7+OPN0S1zPcBKGK5tqbzOZxw5EzHrpGz8t3WKqq60gNactDdLi07pugpsbfZ+zLpt8tC2I2Z9wSJrY83b\/g3hqV0n4cSfno41aCBE6zUrkiF\/nwix6\/KQNfPs1Jtx0w82W+czlwG9jnUX2O4uNNU8+o\/A2NnyHuHuvXEJFwQo4tfMd1ebSv258EjM6XGCNHYiR0Ppyh4KasQ+vg75Xd2TqxDRPritGhXbNid1n9ooXrrPsFVf5zPGchaRsR1Bq640D++nUkXXbq1H6eZ4qzmvsQAls+LvHZ95FHD9aeatP\/08+pEXgLlrkht+fkMd5n8\/BsWnFQEjp8c1RgSsn2vH\/VImrrasBUx3oCnAid3zlweBrscwSSenbtfaCfsuNLVM84Vyklvvz7t2nGQBmsO0exN9aN7RQ9f8zpgC7NfBNHdkHchkxfJyBVK2TZlyqq\/zTM6KOYTJTEoVy2hW83PPYeteL0bF16KBHTZJW9CM52pUQNbJSidr+vv7\/wLb\/gP892fg1q2VwC\/Nga\/aAY25aNxhnyDj8QnAH\/8giaobkYqaKa91J1ayoF9+Grhvgr0g168GejiQuOQQuJgYAWv0L1G0On2+UQsNjcC2L4GHNgCr37JH1bp1TZZ1xTrxlrW1FfqApD9\/+e8t+A57YGSLhcKO1uwsge7AKCwgCJhD3t1jqxf7Ykm7RfYLYs4xOQeBqk9tulnOCxcujMbWVi3Gtk8QhaYZO6u59u0CjZyNzwPjgHFJclCc5+rOVih\/fgOm\/3IbyAEsbhFsDVCEMMfHhpi\/bJGPyxunaVWn\/ThDW0zjHhmCpkp2EcX0BmWikoVG6+oIcvyLC7D8iw9MAhs9s+nj0e72gN0ooz4OTbcKROWs7epz5+yfHgZd6TMRLun+izos2G205sBofWq\/LRDPPWNBm92cJZ++jTs+ejHubPmCgj\/gD60OjJ5F6xHMRBH6qJugvdm3b98eHQ39e+TIkVinYA4FGoQylrPu6CWawCLnL5kzZ2NbEyZMAP0J+5OJDmG2MbQPmwQ63zU+JYdEPFn9eGgBsd79\/Hnu9baLm27hCD1o7SzmifM1B7Ccs2IS0xYjbTUmTNBUkfq7dYx9SkomaT+miNaG6ArSe4MuQ4eWe2tb1KKdHafOaf7ccjYw\/AiPItrFbz\/zGmDgmNjuy+d1mud3wfXxH5lSBE1iTV4MHJV4mE8i6T9veQrLPzd\/0FDTK44eiwGtO3gcvFzxOXPmgP4YHyZnOexSUsqJmEkAfbtj9uzZaNu2rSYTWc4qWM+ZSM5jnqjFQ7mWGNqdx2H6fQWm8zlHcrSbWbSg6TmVibAM4TTr1i1Ex6\/M+43nvv8G5nUfBgw5X\/swEJ0PrroKKO1m6FBEzNTXyGlRz2xrOkqygGihND2yZ9DGSm7pKQ84CLh2hfu2eUreTPU6cTpn9np7gObvFcOAC2Pxa7wB8sClwFOW7Zvr1wA9xHnARSlQo05i+zR1bZdFjX5tmcPehI2UptzUAzeUx1X1vM1N79h7m4DXlkXa6tQ7kkzGkgOdrGXdet60aZNG1EzOiWgugDpuxGwkZxWVlmnkrC1+T5YDLWO5mMmB5LQfx+ECy1XMquvr8UuLrzDrvTUmi5de9IkdjsIJ+x8cu48pirQlmG9r92+P0uPN6SeNW2\/DbgKe2xxfMbol+Ykgd7ON1VpynfkcUWhB02JDwUa8Pm5n3STTuJuA3\/4+sujKxvz2KkcGlRc5Rumkgp2thE6Kmud+a+DKE4BPv445K5I3v8ydZVv4RB9hZDGfNcMRcTuCNh3R0AfspqeA+y3bVHrLiXw0GqRySgqy7KgxoLvUwkdPrWpNCmMt3P90YNSNcXOaz5xD9DLKEDOTc3gUpi1+a9di7Y+vRIXS4v3+vgyjbzEnpM8\/ohofF7jnyNWubzx+B4o+3SY90NLjz8La\/fNMMhijh4kWOL0wWcB\/\/+EKnPrT47H+bCJ5idqZtucCTN9m2eIuHR3JZuX1oUV22X3A4zfK1aRFd9ZKILc5k7UAsQXbXsO46idMv6H7zEW79xIS86l9gL+fmSQJizQnuk7nIWuZ6OrWwG7AyqssnTmlPdWL7lcITCh3zrQmGINdRDW96D\/7nIVD99wf+d9+AzxzO7DkVrk5bCw1qxLoOSD6EyZn7xAGUkOWmFUn50DAS1OjM6tfxcxtT5p6v+bAEZj\/UK\/YdnLLeuTkVSO3q\/luqZvIed9+g8rnHkb+d9+4FcXanv1R2rufqdz1hxyLKzr1j\/5MuwYzs07oLGYk6tqyKscgI9MXR6JBGZ81u52J4h2PmX+YTKhON6cxESJ8Rm1CRWTtXVTYF7f1GCY87pj2v+PZ6V4iwLnOyqYCIl0mYMleXAHcucLc6dBDgecutxFkWI6zhCTDhLlAz1LZkWhhdK98axn++cm\/bOsc9M1XeP75R6XeW2EjhnnM5CytmmAL2jmMiZzFVFdasEgG3zq9pCs+\/wB\/2vKUqTOyes\/5dlKMvBIkZmOjj7\/2Aob3HwXUbAbq3ox4O+tPn2HA+Xdo\/xM5\/Wjb259ujeRqbtpae6NZD\/TcY4sjSLR96eTZbXQOo4aG9ASW1pfE53R+tjF5ZSy9B7j7T\/LtJLDwyzeuRkmRlRe91nc\/UGUJdlN+HjC2KKCxLZwemX\/GJ8EPt0E3AKssvFjaFVg1WSC7dnXvGeC+S9wHlsA8tV5Ls3Zi+rjWj2HIh4Me+r++5a37k1gbmLwYG3MPYG9td+2FqwSTc\/r0QQsfRViiaETGR9uO7lyG4smx7eycgleRe4TZshbeLa3djMX3X4i\/9hF\/xcs6nxi9t7XF4dUNyP8onojpqlVl81KM293sxGYcj9M1GeH29kFVmP6q5X5YoufPIvXSYksLWmVF\/EeAqHyCBJC+meVfz6IoYB8OnoRX3myF024393PpkIjndSCPyI9g+GVAmTlQj5e+rX4PVNfxSpc+b+gDwY4IdcKkoxiLg5ZQNmrz81pUzJ+Imb2Owdbf7mk7hOP364KnjjzHfohOtxb2zUdhVQ47hHmZIOkuqzI5q+wQZnvutLMVzmkoQ9VrrUzb2c2GmhchIckatv3q9tgTK9sV4sJ+4vjXd\/Y8ERcU9LGdfvpWJlnLlc894j5Nm7bQpr9fHLdVrS96dhb0GXcAizaau1hzUhWKH7QQdNNVLHdhPJbQdxD+tRa41SbWN3l5HzbIY8NqFxdtZ1\/duQQlGAgiNusHmOvd90ThEDmAeThntuvW6erX+BJg8kkOZ+Y0Zz77CLjawd3c7aNOQKbklDm2\/wm2JC31cb38AeCO8+KGXfh2AZNzonMwHfWYnNOBOiDKfUterwe+VYa\/lrbChbp3dst6SBEzDUOwiNWddytKf\/OzMLKRW7zftY9NQ9HDM8UA0eLY54TI\/U\/rtawdwJWPxRMuhfekaGaix7q9rZXdVQq8YTlfd1vwklUnLbrzrwA2GJza9DYp\/jfFAXd7Fs8GXl0as8oJnwtuB\/K6K+NsZredPfU3kSh1VmKmhBRJeWDbYWpHzLKWqZuuAJTOAlbbxKKnMVHijM77O4xPtN2u99u2EzBzuVnvbn4Q++aj7vxbcdvuzXFnzUvCEbgmF7E4tG3\/qTn6f3Agk7PEfAhNESbn1KtCRMyN37ZG+3dGoTivFcijVHta1iO392LktK6NCmn75SzyYj3+T8CFd2l1nfLXTjroGNzY9TgzEDbXkW7p3geHDz4fRSWCwCEWKGkRj46l6XeUEGN0zLcsWoOsGGuQlfkn12F0heD+87wYHoFpzy7Jhp2zmOxVlxnPA4dbsA5sEIk3fPm\/luLWD18wzbtzdx+B6fea9ZHQfXsZsfQtZLoGaHwCctaTScDhuuW9cQlA969FD8nd\/mDgx+\/tj1IEY3OK7e1K0CRH03u8cY8uGPnyj0zOMnMvLGWYnFOnCSFB7myFhldGaPGsycu1Yn0s2EfOIZVxntm2L+TyucAd42ODETgzOV3hMFnRAmKibfKyomHRa1YiOfRwhXROqUcba1g3Bo2fm+9uNtrskv9xDrD4ZbM+yOM7\/xbL9vaJE4Dz3OKH+qBXpyxYhO\/EecBXnwIPXu18BmkVxYctWR9GZ9uEaJ6e\/ZtTUPHQ4aY6MiE3E5LT6RpTUEcbiLx3MhHOXEm6zD42vS1xj2yKbGYDmN3H9d+6DcFlHc23K0RNqLzO03iUCd+Z0IS3qaSy0lQ6cxYGHzAQc\/Eh\/ws23WgIzCHYzp7Q4Wjc3H1ovCZFJOKwiL333Zfoukp8d1KzzOnaldH5a998jO3SDQ92MifF7fq7Nri9xwnatZCHtm3WIiCJnl8en2X6ceVkYEBX8YQUpaCsPXwKsMjcBmTjbyf7shC2d54HbF6ZbEvxFiBt0YcwAEqc9\/DOVvhlqTndJxEUfUyKdkESBsppuzcgi1kkK5H0Iy8CU\/7hPBLbXQMax5M3A8\/e6Q4FJWeZME9qHtgFMBmbdwTuP2y4Y18qr\/NMzooGRHef\/ekvIbRYDcSsf4lHnWyImIvMkcKINI3BQEyjmjrYTB4SlplbKjwi6GjwkonzsfbQItgF8HdDmHYFGqrGmYrZWc8i720igPl1gutVqSJoktwtNKhxdEMvBEZcYYr+Lw51AAAgAElEQVTvrKW5FHn4Bn2G7qYcy+8n\/WspbjNsZ5MfhL6zoxf1nZjdjgTSeKWNiPqZauCSB+2BJEc44Xk7kfSjM4FVliA6NB7yXzj1SilSNvZsR9Cntu2OR3ufYSskk7PHFyEMxVVXWhgwdJNh6IsLsMKSxUdb8Op6aVXJAokG47Ah5nkUjWkfwXaZ6GzYA2k5nUVrntodSoDBZZqctFXthaDpg2LrrkgoUs3y2hm7Gjbs9PV45sRjhNDR3VnKWW18tAVwjCUYhMRHiJtuPP3ebrGlRrTYy9cA3UvsF9zZI4G1j8Z36UFfnuT1WFg0F4zzVG\/Ozm9A+73u+U4fInSViOanXYhUN8coau+yB4ESh+tDHseYaHE9\/aXVS11vj5zGzneKQUI4dC9G3eoq1O0fif8tG8aUPlivfdKQalOwq0btCa9W7qrXdrY2bdyE+y65ms+cE50A6ajH5Bws6m4L3sTjgCWvNZ0zC4iZpHN0\/CCLzBhv13ofmBZA8jwm72F94ex7MlDQMxq9q+7p2zD5k2r8Y\/8D48CwvvD0ot\/4\/lrMrYuFGTVWEnmAazsHT65GzVJz2KgZUytxTRfxima9g6p5b4+vA6zneakmaJ2AiHwIT9qW9hKfm6yo28bGT7qSs4HLHgp2Mjq0bre7Y93OdiTmFXOB2w1+D9b+dJwIN6d7wlSPMktRvOyQbfsTSU\/+B\/BoLNV6dJSis2g9WQz5klBEPLtHt7zPOBJo\/Tvg2c3AGhsPcq0NgbOo3jZlumqWk4u8FnuhdEO59oFM72XjBXczOaftDUugYybnBECTrGJ3zhxd8FrWg6ISVb4F5Ba+gpwua+Nanty5GDMPtrlfK7Ka9W1SWgDJSntQFPJIPABy+up4xp+Fv7TmntWdv97\/7ksc16YTVu\/4SLPs81uY44DrjVH5gvHm3+UWz8Xq4aXCHQHR9jZdwyoWJdfwM0CJpG6TKubVEzypzuQqX7zlGdxTa7hsTufMa8tMux0a\/qJ8y04OXHLdx0oNGgucMTV0pGwdhuh2gV6GiJaipK17Nz4KmVc4HMs7ELSo3gHXPY0XlzyvRMZBq\/zsENauna9zJ+jGwu4QduwL81C546MYDE3nd8hpRM4hq01XpERYUY5cypVr+1itZgp0TwHvk1gsiaAnHjcSz+4VT7Irjy5DSevChNX692eBvxp2dXPyq7WrYhRxSkTqkx4Bbl5q7k47q55TBqy0pOlKhwWdMBKxKy62TRw7HjhtckpIym13h2R89CLgjL420rrFnnbDiXQ3aAww4NyUjNdNHNnfe02TKduuUzndwjamb83pthK5B1e5Nt\/+L4+y5eyKUogKqGw5h5mcReez5BhlvLPsNA1u6jYUEzsebV9EZH2R1UwLnZ3zkdu8085NI3mX7c6ih7TphGf6jnJrSfh7kTVMAVbyWsPW2c0anCS6rTptSCyfrd5bGh2HEgJE5sw1YC9lme1s2\/joNpGotDk4YhLQYi+g6uF4PRFYVIYeCibi5VggIaCDrbTyTWCwZPIzLY3mPsBxPYGdPwLr3zWnTrVKSuUHdQNGHh1\/Rq2fg1esAyreqQXtRAkfMgreKUb+wtOYnIOdCv62zuTsL57UmnA7uzEXyGlw7UxLE9llAEa1jziLCR+RZXzypcAJlwAUtMGa81Vf4PUFkc773tsIvLY8dm5KxK6foTZ1aucZKhVC0EZ061UpWlDog+VPBUfi9p4nxNUSZa6KenpPGQS8vspcRzWCJullPMENH06uk8hDgcvefA5zPjIcoFq8s4XXhWievPIscO\/F8di7Re5ychDzIHfYihJR3vI8cMfyeMm0nNYnAsd2t48yplvCdV9G6hOBa3+3lhupTtRjHqpH3bc7kdvhZeC\/u6Px3120GAr0FLxUyOQsB2c4SjE5+6sHt1ytpt6avJfzGgoxtX8eClruLfbItopoXcx1QqK7uHQn1\/hccAcw7KKEB+mUHJ5IekGvP6IBjZrzid15s7Fzq\/VMxKx\/8dttb9taz9TwtScBm56OH1\/IrihJKcApcYHeQJsC4PrVvmz\/inRr9c6OcwAjYn7hn8A8S35FFT+KpJTirRCRJBEqRRzTz+dlCdZbT\/alSQYtrashzGrzH7fjwOr+TM5+gZyKdlQm51Tg46UPaWIm6+TtAdGrVI4esCIBrGd8p18NHDowsp1tfCjAAS3kST40rvHVT2CNTZARY\/NE2BM7HI0ee+5n+6FhJVstdnjLyJUrq+MZ\/eyaf0aukuhP3DbrU3OABybGjzIkV5Q8w08E+P4m4G\/291Zx3cqIzhN8ZLazB3UHVlxp6MDOl4GJOUEtBFstak3fB\/zmm434ce1IJudgIfe3dSZnf\/CkM1pK\/6jf67Vt1eIFe+MZwF\/jd3PthbrrAuB5wyexvjBaz5kDcJByuhMtEnjQvh3x\/FHxMbhLrzfc2aTQfIYwpXbBVqyEHuc5bLc1PPpG4I9\/9UfJqW7FLr60LkcS59FxqSAF29mmLFM2sda18+KQRjpLtbrC3J\/q6zx7ayvmrZ2Wl0G\/o0kBBXbVY+9f7Y7D1tzpTsokrIWYbR1t7AYmIiAK07l2YXyEsAAXzPvrXsafXn9KCn4R2Qodw069OtoeXckiC9r4XFIB3LEi9hMtatj5FhHsriipdtXKiqzbtTiPBCnczn7tZDTW9I72TIFxpvdrynttl7\/YY79SE4YLBYIAk3MgsAbbqMpKS7m39tiCaPCEv\/fsi8m9LdvILqqiEJa6c4bnbD6iLUXdcrJuZweYHMA4RP2u86ovPsSLX22NJruwwiAiaKetbapfsk8hVvaLRCajR0TowrCJmWzhyZxJu1jTwrCtltjZ2kfj6U8Cs\/7H0HZP\/zOAKwTRzoJdrrj1BBFQeZ2nIbPlnArL+fM6VH2Vj+K9m6IrJTjZNIXlmMM5NjY2JtGaTVUiRbIcKhdoBegecOnxZ9smQyfruPH73yFnn62mBq3ETFafMKCDSAzRlRvdaglDxCyDzLT4\/\/n1p7DcEq703PaHYV6vP0ZLWtNJtuvyIT7tbo5BfOL+B+OJPmdH61ijht1yNnDpEAFgTleUVHQUsw7xjTXA5AHuc33cLUDf4Sbnsbjt7P+d6ZvmZkMd5n8\/BsU\/29ybTWIr3V1gLhEUAkzOQSEbYLupUJrmmLCmSiM54UtP245k7XkM1Rc4OVus1YqDemhpE4UPndm9PhRo9SlyDzE7YTV+UYCGtZFtWs9JA0iGdY8BFVeZu6X7yKsqzGEQQ7TNePari\/DY9jdMMhudvUSW8IAL58ZZ3+Pye+PeQ0\/W2hHl3LVLoKFdC7v7T8Brz8erK0Q4Jfxqa05jLwN\/O929if0KteA0db\/dEx1X3GT+aNwyBI3vx1IOEjGP\/m\/kQ9T0EGanXAb0Hub5PXUXkEsEjUAq1vkgx8CWs1+WM0VvWv2Qdt\/2q5y9sXfjV5re7tztIrzRrAeO\/\/lZnPRT\/PUXY\/hI2go9Yb8uGH5AV3uv36AsZ\/1c2bBdbA1tee77b0RzGxe9kYsFP12KnEPWIGdfQ0QwGrTlnNmTZ7bFao9O\/v07Avu0i7\/PHDKrkJJk6HmdSXbr9rZ1a5u2qcveiyfoYft1wZIjI8kPrPekbUNK6mAtux+403o43fRL4\/1vjx+GQS5EntrW5+qtYxzjVVvzcWt9WLezG+pQ+60luUomYOQJ0MwszOQcMr0WFsZCLdbUiHPt+qY0G8\/Sutx8jNl9Pqr27Kmhk9Pi6yhK+Q11qKmJ355b0KknxvU\/PlqOFvXyXqcgv7EQdEn\/kRciMWvr5njY1tYDc\/xrHfCvtUCzXwF9TgCOiPWjWe6CLVF9YXu59QF485\/3a1vaM3sdEyVnW7U7EbMe8IP+Nj4kA8lKO\/STbc60x90MzLWE9bzgTmCYOC52uqalyPGIAqzoyS6sW9u6k5eV1En+4n0Ksep\/Z9DWoCRaQowpLiOUDfJR2BM48uRoQpB04ZZwv\/o7KHDgEu36GO8007to2s7OhN2FhIFUo6J+VYreAdqJMj7WjFe+rfNpgiajLGciZiMhW\/+vY5yI0vSL9lE92ZzxVbTtj7Gd\/4KcwpddVZr37Td49\/F78UmLPZD\/3Tfa2S5tIa\/dPy9WV5Bb1thwYesGVLYsAy002qNbFa69OxdY3q4Qnb\/+j1bI8bzZ2IzFKikfD4w9uA54fSVQtTDe6pWV8cRLgKdvN5cO4NqUrDhu5UTXr34+OZIPUrS1rW9Tiwhat7xdr1WJhKK5cP8EccASa\/lMOFel8T46Q4s+Jkpo0vh5RzSsi11zi25nZ8LY3Salgr+nD10K8kPvzIPrgfnxOXKEo9LChebVY+Dv3kP5LL7nnHbV2xGx6Ocy5EwTgvINU55d09dZa+DeEz9G54XnYmtODd7b8\/dY1OEQd4vSBSEi6nteeB6dmghxfP9haMghZ6y9kP\/d19hvwwAc9MM25DdujZ5h1+YWYEDL1Rox3\/jDlTj9p0W+6IEWNvpYcDxvtvZksZjn7TEdYz6Z4Ys8cY0oYOE0XxK7JkXyG61n0da2HlHJjqC71ZTh6fWxxBzCa1V2aBtIy1UhCmDrOobP6zDwtUWo+u\/\/xYpa5uf0H6Zj2p4VkTjX5P\/BT9oQ0I+BXvv6Ezz32buo2\/W1+JpmUyAjzQG1KTynSGg9Ah\/ljs+7\/XIOQpI2zTZ17Bc5k4VMW4\/alknL+tiW9D51yNm3hvzbpRM5JIoJEfUNr6xGmy93w6Kfzsc9u10obKpLw7u45\/sLNbKuzc3H2N3n4+ofr8OvWnyJov\/bgiv\/MADf\/mo3\/LH2HRT839fArlb4JLctpvz6OtDWO\/3RLe6rfrwB5\/33fuCQo4G3X8Da\/dtrFrPt0xSGs6GmNxo\/7mFKs2frYOMFkIN6Ax8I8icrQh5O1vOEB4HbDfGIrUQ7cEM5qqyRyQTpDIXXqpww1reAd\/1fxKK2exS3JIUZpyxX+mqvSv7mhJfpnE1l9euGL321DSu\/+MD2umGymDT++2A0vn90jKhpvd77Y+QeGTNSDrhqPV5cuZBTRiYLdjL1EyHn82+fpSnt+9wGvNZQD8rTa3TmSUYeqktbkvkt9kLvVu1Qsk8HbJw5HTPzzgJa\/Vue4F22takfItiuDW9h0o83aYR78VHHYp\/cz7Dt+x7Ridu4ay8TgVrHRm0c9\/MyjeyFOY6JHGhLsDEX0\/4zBzN+Mz0OnrgzPBGAejYe\/dxZd0rasQ3Yqw3wq1\/bb38rQsz6sK3Wc3mvEVpyDxkP7CEvzsfKLz6MQ9B4BciT9SzSBengs4+Aq21CYobM2U7mfRSF6Gz8ohANlKe56eaAKQqYTKNcRhoBL8F6pBtNoiCnjEwCPL+qeiXn0y4aj39POdGv7iPt0HbLZ52Q90N3rL5UkBTh8zrUfbILBTcfErXKc9pvQU6hwEK0SqaT9BxLbuGR\/txzJo9VIlfhnWbqu+4wNL5daosX1V2zsyR29q2XJEI9dhzw+7aR6yjkuGPNICWjhbKbgOEOeZ5l2khxmQVbX8O4zU9EezV6bjttbesVbni\/ClPfXhkn9S9LJ0U\/tDxbz3Yk\/dhMYKX5zrVWVLFY3X99axlu\/mB9bJSGj9sObYAPbwl+EuiWo7En2SQpwUvnfw803mc\/excT33jW38abduf0GPReG6cdyMZpS3lb2ytwfpdPOTkTEe\/aC41vl2jWpDEKFoUBJKvG7ql7qw4F1zfldm0qpOU83uNL5B6xxBGaX4Zfb\/69HTlTUoWdrTQLXZOt6f9Ojefvvhmdet2EygPMslmz9ljbGPzzCizfeWzsx0TIZ04F6K7p+n8AS+9JXN2Kb7FarWf93vP4ucDcNTFY7Kzgmp1fodPKm+Pw0y3opK1nY8t2kcYUsaCdMk7dPQa4MPGcGVLz94G6V3Dh6\/bvL32cddrj97jn0JOlsplJdWpTSP9AqP7633jvux344LtIXkbaxSvapxAtmv1Ky6qmP1t3xW6U6D+T+aB4dPsWnPPqP6REpfFH+9sB9GqzF16ryUHjZwcB9W0j6+iuvbSjRFqzyA9jfocFKF4+RjMa6Kn77V4o6z\/MPiBSUwdEzB8tugv9PzgQCytf5G1tKQ0FVMhXcv7+d8DuBmcS2qL+9hsUf7oVU6vXRzyrm85s383tgrt2+xP+1ay752AbVuspCg2dnbR7C7k94oNJyJBzTn41cnsvjkfakhnKWsCYjEH\/nRMxay\/PyXUovslyTzRZHetb3xR0RHFnHbIm7qx5KYrI+QV\/wF09T5La2tYr2SXfoK1a2rL1xXrWO1OYoOOc6ZpuD8w7DxhTlOykFNcnEiRHptNf9hbWk4hqSucSlLQuNBG1TqoPbqtG1Ze1cY5RVK\/i8D+i3e57muppDqzL6jUn1a0tN8sfm0nAopPqiAO6ovNvW+OLH3dilctZMtWZ02MYftv813ExG\/TUjno+Z6sI0TC\/WxZEcrVbn33zUTf1Cez\/3TcY+s2H0aNIWqPJefaa6vUo+nSbVqvw7QK2nCV0HGgRz+Q85gpsP7bCZFHmffc11nwXuWerO0PR\/V4KvkEKJ1I2PtH7zM2LNWJec7V8onC9HWuACf3n1N7o4+rR\/\/B6DHqhPNqtlZybPTk5sp1OZ8qNufEBQUSoC5yLRIROX69koYleniuGGSwRu+QLThonAqaHtrq7HgN07gP8piW0jD+qBscQjJd8GIg0jI+ew1lma1uvJ2qHftewbjQm9z0I153q4+tlp8\/is4BJD\/vYkX9NifChuTvvlALHXaxkJJBOl5pMJ4rVJVKe12uEMIgSkfGmj4Az7hAPKrqG0rXQV5cC91hiGNDaMHAUMDLe30VrUb9G2hRPYWPt5xg58Som5zDMIS\/3nOnMud+AIej74Q\/o3vCmJr5tbF3R4PbNx+ajp2HLQaN9efmvfwpY+SZweAHQ7UDztrhd7OakMd\/ZCv3eaI51v+uO3K6WHMgCAqeXZ8rJQGlXwUcIvRDkXX3jaWKx6MXqNRgoOlN5a9gr7h1W3GSygHRyHvcAUG4I5+y2RS3attUIem0ZVp9XKB+3XHYA1hzaVI\/0SDsaA0fLthJ4OVFii8ZPuqJywEj\/MWkajZ0u6NdEUKWtO6Dd7r\/DKQd0w\/21L+Pu2o2B45DODpxImeQiYi6ZFfnb+tC6MuoYYPqIJoJ9dCawyuL\/kMDxlsyV2XRi5tZ3RgUhocHKRgjTHcJoK2TK5g0Y8\/4WN6xiv0+T17Cd9WQnOL0wI9v1RG5ODma9ZzjgdBkp1Wv\/VplmkQ\/oChzTxUPCCv0LlqxfilJFf2eYNSw\/USIlRXqjoCQyXtuivqzn2FTmmM\/HYc35Ph8vuO2IkF4pcttv907rB1fctv\/OVrhh\/zJcPiB2xulVZ07l7YjZlaB21WPK2g\/x2I\/OfiXWvqld15zp1kpNu2l0LKWd4dLOmh6tsGU9zjmmEQV59fjmpx+0n7drsSfa\/HoPdNljX82565fGBsesa\/pHCN1Goeh3lPLU7iFCph1CESmb\/HNEWeioUoLrLZOzn7M+RW2R0qze2kTSk19\/AWXvvW4vBU2SkyYAJwlmWopkp4Vh2efv46It8XG6dRHsFgkZcqe6lf3KAndYSRFcoenGzjEskchfduRQfuAkjDo8AEKieOsyHvZE1Cf\/JaWYi7aW2+8oRs34QYHJYd0JoY6WHz1Ws5adCOr5LcCfdIOQ\/Er2ex85B74ZOx9u8k5u\/OQQaHd4m64\/nntM5MycCHrhC8DK2q8iMe3JidTw0Lt7ZrseOCjnII0s3\/8MOPZGexg0n5HzIrtgehAcUWnCmB6j05iMsxiRsn6+LCJm0zFgAGlPmZwDewWCa1hX2rVPVKBt27b4qfEX\/PuHbzVCKvolJ7IQVVZEzjDoIeuvdFTEOkjzeag1K9WA9Q9EgTqtbQ+Qw5Hbc\/7mJ1G+9dW4YkzMbsgl\/vux1YtBTj76c277XtrZ3MUVwJ0rYu3+\/UzgcpskYMbe7T609C3zxCW1qSkTq1uvmkLv7gVbqzFuc8z5Ufsw7VyGojz\/P1KEeaHptlnT\/XU7zCk2esV6+y1dSgxDcaEphj6lCPX6dNgXWDVZTLBOBGnsh8iZHOcKXIhaVjZTMCdLpZvPAv4ytOmHtMZuqQRuj\/dtSdRi1rtjcpbVVojKqaw0v1JG6nFrr3xrOXbLzUVp646OW1MhUp+Sotg5htVtbRW3INumhLSMXNRmXmMhPhoeCbjh+6MnLpGxpCeUA4PG+i6CsUGawx0eWWyyIKd0HoDpB9vfx09GIFFeaLoaV\/t2ASrWxSdicOpLs1oF+c2J1DZ+CJx5p7ukXlKxUruVbwHk5yDz6JY0\/d2zfeRMuFXLSE1HK3sH0NAIdLg0vpc4eZ3Cyp53G3CiQxQ7iUGovM7T8DLuzFlCZ1BZaX6RswxOXMZfBKzbofqdZ88pIQ1iia5ZHVM7CWsu9d9yFKJBu0x2qRvJcczOs9YHaEvuq8X6NjFP+CB3fkRHCV2+OwZvPX+c55HI3Oxwsjx1gkzkdgi1W\/MFUGoJl+B5EAaSJmu7sSkWkjVTlN5u3Jidzpd9ineu8jrP5LxunXKX05mcE1lGwlHHaunSuSARtDWN5MWDgdtHyct80saHtIQBxmfqbyZh2nEpImjqOMX3o8mZjvJgG89d9aMCeeTkSgq3sy0Z2GRa8mLpGtubVwU8tIFMKeCJicCWbR4dNG2E0+5GL\/Zm8cuM01pm6KHAXaMNFrfPjl92MjE5J6KtNNdRWWlMzmmePEl2b3UMozPiZLa2SRyRVdf4aWdUFp8b2FUiIQx23t0+hwDVruXcWo+P\/3CTSQx9JyJJFcVVd0ukYawg2vIddljkiqRT1EC\/ZXZsTz+eoA+q7sWgiIUUs+GbXcBtz\/tH1sKPkRQRM41f5XWeLWcFLeeUvsTcme8I2DmGWbe26TzSy2Je8U6tZkkaH4rB3Xh\/Cq1nOwva59zbtz4PTKp7EDn7vxcdrr4L4bvCAFg\/qBo\/7IuGzTGvvUQt4iBkdSXljUuABwQHwpb768b7yOSoRg+lz\/3sG+C9f9sTuP5xIjpP1xxsyV\/hiyZHW13YK\/8B9PMzik6kYSbnlM+w5DtUXWnJI8AtpAsBuzvP1q1tt4AkIvmt4SspwtuQr8fh2ctTPFqRBT3kfODP9yYtiH5nttmp5nzZQVnN095ZZY4RYAnOc8co4KLBSQ8r8Qb++z0wbWjsqpseeW\/m80D9Z5EbJqQPCoMpcx3Ow51ia0ARJ0cxW2K+6H7guPGJj9+hpurrPDuEtWsXyMTgRhkBOwSsjmF0Fado915xgRq8xsy2C2G5+vyC1G5v08BFHt1JXrHSnaTW5lSaItoZs335OeuEiTTeGqBlZwuFtTy2IN4K9QsAJ13pYTJrXgc+qwHadQbad7UPRGPnjzB9KXDEEL8kjmuHyTkwaINrWHWlBYcMt5wKBKxnmPqWrHVr2ys5k+wUCY6sPeNzYNUsLTlGyh8RQSdx\/qxFVHu4Hs2KyiMx8Zset3vGiY77ii2rcEutIbJekxOYjLd1on1K1fNy59yuQbKQB48Fdn4NPGmTR9NqRRPJuqV81cNs6hb80ruBxbPjpZj+HHCEftlZatSeC6m+zrPlrJjlzA5hnt\/RUFYQRQyjbWhjEIpEtrZpsNa2G94agPJ+pYFlZrIFWLS9TdusZJUl8FA0NfLOzi02X5\/6aPCkBFpzrlLxWj3KPjY7nFEijTwUpOdDRxc3WWIm\/K3Z3khP720EKAJOKp4U3IGnYTA5p0KZPvehstKYnH2eDGlqbuiLC7Diiw+ivV\/VqRjXHjIIXjJV2Yku8i4m57DaWa0cA0gEAoVoS\/PUq4BR8hdtjXd+iZiN16du7HocJh30v6wJPj7aNa135iGnzUexVne2wvy8STinn48deW1K9LEjSkRCmH9WC2x9E\/hoM3BwX+CATu7JSuy2n73K6WStpzBpisrrPEHIljNbzn69etyOBwREBErXqsbc0kqLR6w\/a6YkdqfVeq7d8OZgTO5S5G9qSdnxUqayDY+bS0+c704WTTXmrgHGk7Hcsh7NhgZ7fYqIeewTZuucxJj6+3GYdozPiUVk8aNyIqeuIDKEUT9PzwGecjgHMeZb18+f31zr7HDmwdHMCyxOZZmc\/UIyhe2orDS2nFM4UQLuSuQYlrezl2lru\/hggAja6xMq5zASPkEHMWPmLkoUkXtMRRQKv69P0UcRHSvkFpUjZ9+awPrxqkutPBHm\/ZaEO0k62LnKsew+YO1jsRwDhx8L9D\/D3vGLiJqs7x++A15oinfec0AkXzttp6f4UXmdZ8uZ7zmn+HXh7owI2F2r8mNrm\/qxkn\/jF4UYm1uGBwQ5BgLXjF2AEgeCMaUaFFjNfjqC6X1Zz7QJl6CuaUljLtpuPncWcNpk6SaysSCTs4JaV11pCkLOItsgIHIMm35vgWlr22tAEr0r4VWgqnFIy9UqEspDiE9rbOmc\/Grk9jZnn\/LLESza17Z6rQ\/jmTalgaR0kGl77M6Z55nTRaZNvhB3rPo6z2fOip05h\/hdYNESQGDGu5W49t3Vpi1UOt80em0nurVNjVodz7CzFcb9OAn3potvJAn6nlWG3MeA5qFtJM1rugzANV38yT6lb52Hzmq2O2em3YY0p65NYKqnvAqTc8ohT75DlZXGZ87J6z9MLYi2tmkbtXiC2fkokTvPNM7QnT07WdDXrwZ6lMC0nU1eq4LrU5X9yrT868k+0b5o29xyf3rSQf1xY9djk+0isfpEzIuuA1aUm+t7cKRLrOPMqaXyOk9aYMtZMcuZyTlzFg99JNawmxM7Ho3Xlgz1ZWub+rir5iVMeONZE3AX\/TIJt41IntwS1obNGXTdWfNQ8OwYU7NWq9nPrebSWcDqtwHRtrlfHwCeMSJsKhdEAn4Yn+GXAWVmb3XPbWdRBSZnBZWtstKYnBWccC4ii6xbSvk4\/aEYebjkD6QAABm3SURBVCaztU3dxzmH7SjARyePS\/29ZyMWNgRd0nINqpo3efcKHMHoypkfVrPunS26ojW18wBMO9ifbXNPM5YwESWn8Dl5iCeZFC2s8jrPlrOC3tpMzoquFC5iW8lzSucBmHGdmRwS3dqmrkUfAOc0jMP8U9J4d1fHZFhOHDoFv63V0hgWDl2MrS2ro7\/3K462lnZyFrQt9NwezyOnMyVMjjx+9eF5pqYwnaJn2RSswOQcMqUVFhZqEtXUxO4pWkVUWWlMziGbcD6JIyLPA1+O5HrWn0S9tvX6AzeUo+pL83vhlxWaFAw2FvRfhlXj9v3+aWrar+tTH30BdKTMiQFf0ZLCRb8fTJmjrA9bzFIQigqpvM7TeDLqzJmIWSdl478ziZwTnqlcMfQIxG091\/VCwysjonInGmtbb0B0tar0952x\/Jhz04oNWa8zZlRhfl2JSY66PfZExzP+HP2ZXxat0eHMep7tVx\/SgDqlc2RiloaRyTkpqIKtLCJjO4JW\/YsqWCS59XQhILKeKSY2XX\/Sn2S2tqkNO+9wiraVrufCecC9lUDxz1VYs9NM0MvadcCw407XRPPDaiZiHns\/sEbgBEZ9pHQnwSmW9YV3AsfHPkzSpRuV+1V9nc8Yy5nJWeXXiGXXEXCznpPd2qZ+Ojw5F1tzzEEsfj55VlqUYAzPSQJc+8NUTPnxOpMsY\/sPw9pDi+BH0BGTE5jl6tQFBX1wZ88TveGgx5be9DTwy89A3+HadTDX59GZwCPT4ovpKRfTEO7SVWbFCjA5h1hhbpbzwoUL0a7pKpX+d4iHw6JlAQJC6\/nxGHEmu7WtWc9b61G62Xwl52\/dhuCyjqlNuWS9z0yyUa7k2t1L4pIorL1sHopKBGeyHuaEsb+cg9cgt1ss7zVtZ3u6OuW0HU0yEcmOvxUoPNQcMGTh9PgrUvoYzrwGOMtyfcrD+LgosH379igM9O+RI0dinYKOvzSIjLGcrRNT5szZWGfChAmgP2F\/2CEs7BpKXj4367nxkeT7GPBwNdbtEQuHSS2mOoY0RUEzZuAiYl5waT0q3p6HeXdZcjQnef5q9M4WRQLztGW+eDYw\/wp5JZDs9HxRJ64TRHYpeekyquScOXNAf4wPk3OIVOxEzCSmvt0xe\/ZstG3bVpOcLGcVrGcm5xBNtIBEsYvq1bgjci7sx9Y2kVWHJeaQmH5neXKCZ\/piYMYT5hI0LgoGUla9GHnffoOPFt1lLkBbvRS60uNjJGbyzrbGz5ZyAnPyqPYoj6n4foXArEoOx5kMhoa6ZC3r1vOmTZs0omZy9gncZJtxI2YjOauoNCbnZGeIGvWtUcPIKUxzDgPgx9Y2tVNyXy3Wt6FEybEnFdaz9ZyZer\/seOCik+vRcUVsu73Vj99jx0O3mhV22CDg2hXSSjQRM20VHlKJ3K6xWObUkO2YiZC1AguAVRX2lu\/QC4FOf4gQ7Hf1wKxT3OVLcifAvQMuwWfOIZoDMsTM5BwihbEotggIM0qtLQOlfaQnWa9takOznv\/xMHLavhOVQ8qKTEJvdufMa66v1yxm2jXQH02WqqXAGxZrmYhNIvmD6COg2alXm6S\/pfvxuKTDUZEz7tUPAVtW25OwddxOchCxv7YMuPtCcy12+Epi9niryuTsDa\/ASssSM5NzYCrghn1GYOrbK3HD+1WmVvWrVX5sbVPDJ9xdj+cPMDuHBWU9W9NAUv90zrzmauDB\/1RipiE7F\/0ueq3pqngHMQ0UIjqyVrv1B7oeAzRrHv1ZHDHbbWfveygwZZA3zSVCsETWnEnKG85JlmZyThJAv6rrkcGs7YkihamuNL8w43bCj4DVOUzf3vZra5sIs\/COTcjt9bQJjCCuVtmdMxcfYd7OJkGuO2QwruxUFJHJzTPaosZ5u41F2e7mbE7WxBZUpfK5h1H06Tb5SUCkfPpk4NBBTLTyqKWtpOrrfMZ6azvNCNWVlrbZzh2nHAHR9jY5hjVUjfNla5sGNHkR8Pdfm7d7yzuPw6iD\/QtMItpi\/vMgYNJpNtvZgy3e2iToynnAnDJHHSzYbTTG7D7fVEbknT39tXWYsjkWTzuuUd3DmqzdP90FfP0FwHePUz7\/k+lQ9XWeyZlTRiYz\/7luChCo2BbxYDY+jR\/2ReXxw0DZqvx4cifWgkJZGh+\/omXZ3me+DRCNzXFbnazoFxcD5fHkTZmsKKOV6RHkaRZ6gutb1VSZSdiPKZX2Npic064C7wKorDT21vau70yoMeylCiz7\/H3TUPo3DsDq4f6kNaT7xgOrb9ISQegPZcaa7kPaxIKJEecz\/dHPmakvo3c2\/f7hI07HGe16yKlM96b+ok5rv6C8Kc1kU+38hjo0L5qL2v0ao+0RMc9b90xsOzuR82M56bhUmhFQeZ0n6NhyZss5za8Qdy+LQNz1KnKaOnqWbzmZRYFJkrWeB98IrHzTPELyNCdiFnpni7azJQCi9I9Vb5s\/AEp+dQkePPq3ptrla5\/FqA\/eiDiOXTqfrWQJbFUtwuSsoOZUVhpbzgpOOJ9EFp0\/U9PJEqgunsh6Prd9L8zrFcuM5WUoonPm8vHA2OIEtrMdOhb1U9w\/\/g63djUrdx9g4Ggvw+CyiiKg8jrPlrOCMVeZnBVdKXwSWxQ9jJr26\/rTZctrMef75M+enc6ZRR8ZFxf2xa09hnlGyRpkhBrIz6tHbnE5tu6KbdF7jp3tWRKuEDYEmJzDphEJeVRXmsQQuUgGI1A4qxbbupoJ1M\/gIc0eSC6sp5Awm+4z03nzJVuewd21G6MaSoY4RfeZC4csNmXdSqb9DJ5GGT801dd5PnNW7Mw5498oHqArAnRfeOb6eO\/qYft1wZIjz3Gt71Yg2ZzP1\/wTuPZJcy960JRk27bKnnOW+ScFRZXYtq85PKenpBZu4PDvlUGAyVkZVcUEVV1pCkLOIvuIgL5lLAys0W8cKIFFsk\/phnKs\/bLG1IyM85loO3vyScCs0wDazqZ2rdvNieZpnvo4cN2SmIhklX9cbL6vPXDfjlh2VHKpJpPFkuunBwHV13m2nNlyTs+bw70mhcAfbwMWvwLk9qtAzv7mK1Z+nD8TyXZ8wUx0g1p1w\/NFZ9rK7XTOTJVmvmsO0UnbzeRslsjHRFxfovvMXnM0J6URrhw2BJicw6YRCXlUVho7hEkoOAuKkGc15UOmK0nNhppjY\/uV+lEUIMRui9jtnFnkBPbAYadgTN7hCWnrb88AVz4Wq+op21RCPXIl1RBQeZ0nrNlyVsxyZnJWbYkITt5hNwHPbRYTtF\/Xq0R3q6\/pMgDXdIkFPxEltKBR6+fMfm9niyx0a7apu3qehPML\/hAc+Nxy6BFgcg69iuIFVFlpTM4KTriARI5az\/\/bMs7tPx85bT409eRX8ormS8zb29RJQYtWmNtrBPIaC0ABQIwRwOj3bufMiW5nU9ul1wOr34oNlcKOUvxs\/fHTcz0g1XGzKUBA5XWeLWe+55yCV4S7CBIB2tqlLV7R9nYyAUSMMtvdrdbK7GyFhrcHoLGuV7QKOWZRFDCRxUyFkvGeNn6QaAtYfjVye5vjjvtx5h6kzrjt1CDA5JwanH3tRXWl+QoGN6Y8Avp1IhFR+bW9vXZrPUo3m8+2hcDtbIW81pHfGL2y9bKXduyH2d2GJIR53Lm2IEezX9fJEhKQK4UKAdXXeT5zVuzMOVSzn4UJBQLGM1jrFq9fzmE0UCLoAYteQk6nFxIad1neEbjvsOEJ1aVK1nzQ1o8RDjaSMLQZWZHJWUG1qq40BSFnkQNG4Iw7gEUboZ29WlM\/JrONbBWbPgRW19Vj3LI3kdt9ufSo7jt0OMryj5Auby1o3c7WtvGLyk1ZtKyOagl3xhUzAgHV13m2nBWznNkhLCPWjUAGQVeriMRy+z6KnHb\/MvXhl3OYsVE6U6bn1fpPcE\/tRtDZtPEhS3b4AYfgosK+yG\/RKuExi7yz2WpOGM6sqcjkrKCqVVYak7OCEy5FIhtJzHq16O\/dhuAvHfulRBIi6bwWeyVFyLqgdvenP+5tzj3NVnNKVKtUJyqv8wQ0W85sOSv1wrGwzgjo27+Z4MUssphp9CVDqrFuj5iHNl+d4rdChACTs4LzQmWlseWs4IRLscg6qal8\/1eUo5lgLD8PmPXjTSZPcD\/P1FOsKu4uQARUXucz3nIuLCxETY05eD8NWmWlMTkH+DZnUNMaQV8dH9rztLY9sLD36aEdKcm96l\/AeHNGTE3e8vFAs8JqlFWz1RxaBYZIMJXX+YwmZyJmejKNnEM091mUkCMw5O\/A8m\/ig3SE1dIkYl76OvDnBWZgKajJtFOA4iPis1rde+hwjEvCCzzkKmTxkkCAyTkJ8IKqqlvMmWg5B4UZt5uZCJAH97rW8aE9wxZFi4j52c3AxRXxxEwxuosPjkQc67giFgiF7zVn5pz1a1RMzn4h6VM7RkJmcvYJVG5GWQR0BzGr9zYNKEwEfV8lcMG8eGJeczVAlrMoFOiTR56DE\/broqxuWPBgEWByDhZfz617IeeFCxeiXZO3tv635w65AiMQcgR0B7GwErTI+YsIWSdmgtca35ut5pBPujSJt3379mjP9O+RI0dinYI5FGgQGXWVymopu1nOxvkzYcIE0J+wP+wQFnYNhVO+s+8GHnk93kEs3deQ4iJ\/IWIpG4lZlAt6SucSTD94YDjBZqnShsCcOXNAf4wPk3Pa1BHpWETEbuQ8e\/ZstG3bVqtPlrMK1jOTc5onmsLdD5gFVO2ID+\/pZ\/xtL\/DYBRgxEjO1N\/PdSsx8d3W06XR\/UHgZI5dNLQJkLevW86ZNmzSiZnJOrQ7ietO9s0ViWD22VT6LYHJO80RTuHsiwzH3AWt\/2ojcwyjPZOxJJltUIpDospDlrD9kMevOX\/rPRFZzWL3NE8GB6wSHgMrrPKGSUdvaVjW7Wc4qflExOQf3MmdDy9EAJUc\/hJwD3jUNOZUOYvevBs4vNyO+ZkrEK9v4lG6Ya4rZzVZzNsxSf8bI5OwPjoG0konkHAhQ3GjWIUA5oK0RxAiEysMmoSgv8SQVMkCKwnKeNwC4r8xcu2KbOeCIJl+\/caBteH4YATcEmJzdEErj75mc0wg+dx16BAom1+PjP8TuDZPAjZ8dhMqi0XEWrF+DEREzbWfX3mbuQbSdfVybTni27yi\/ROF2MhwBJmcFFay60hSEnEUOIQJElONXVmNNi1g4TBKzoWocrjyqANef5q\/Qsg5gRMwUotOYgpK3s\/3VRTa0pvo6n9FnznYTUHWlZcOLxWNMHQKDXijHmh2xGPSNOwo0gj6tD7DoEn\/kEBEztUwOYKP7m\/uwemfTb9kJzB89ZFMrqq\/zTM6cMjKb3lceqwABa4AP3Xomkj61D\/CPJAla5JltR8wiWVLtSc6TJDMQYHJWUI8qK429tRWccAqITPeIyWI1Pr88Pkv778qrgIHdEhsEEfOsp4C5a8z1LxwI3D0mvs3mS642\/ZC3sxPDnWupnX2Q9MeWM1vO\/B4zAhoCVmL8ZekkYGfEc1t0zckNNiLmRRuBKx8zl7RzACvdUG7K08zE7IYw\/94JAZWNMCZnBWOusuXMC1JQCMytewUXvL5EaD3TD8mjmohV9rljOXDJg+7ETCUueeMZ3F2z0VSYr03JIs3lRAgwOSs4L1RWGpOzghNOIZGt1jM5htHZMz0UIIQsaJmHtrHHz40nZmtoTiohus98\/2HDMTbvCJmuuAwjIERA5XWeLWcFLWd+DxmBIBEYW70YD26rtrWeV1wFDHI5f5bJMqV3ILrPzNvZQWo4e9pmclZQ16orTUHIWWSFEHCynmkYTtvb0xcDM56Qs5iplHU7m4h5Xq8RHAVMofkSVlFVX+fZIUwxh7CwvggsV+YgINpm1j239VGKCNqLxUztcFKLzJkzYRwJk3MYteIik8pK4zNnBSecgiLHWc9vDUDj26XRkRgzSJFX9oo34xNZWPMyG2HgKGAKTgrFRFZ5neczZwXPnJmcFVshFBVXxnqmoZ1+JLDpI4AI2vg4ETOVs7ZP29mV\/cqQ3yLYpBuKqoPFTgABJucEQEt3FZWVxuSc7tmTPf1brecTmg3FkseOdgXAjZhF29l\/7VSEWYcMdm2bCzACsgiovM6z5cyWs+w853JZiIAolOY5DeOwYLF9ykY3YmarOQsnUpqGzOScJuCT6VZ1pSUzdq7LCHhBoMOKm0xRuyiXMgUHKZho3somUp46HBhb5Nw6O4F5QZ\/LJoOA6us8e2uzt3Yy85\/rZgEC1u1ta4CQqnciAUpkntINczkVpAxQXCZpBJick4Yw9Q2orrTUI8Y9ZjMCZ726CIu2v2GCIJEUjiInsw8HT2InsGyeXAGOXfV1ni1nxSxndggL8G3mpm0RsFq8VNBL7GvRdra+Rc6wMwJBIMDkHASqAbepstKYnAOeHNy8LQLW7W0qOOKAbvhbt+McrV+7O818dYonW5AIqLzOEy5sObPlHOT7wW1nEAIi61cf3pIjz0G337WJI2mqc+P7a0EZr4xPItviGQQlDyUFCDA5pwBkv7tQWWlsOfs9G7g9LwiIrlcZ6xtjY9MZ88x3V5u8vansue17afGz+WEEgkRA5XU+Iy3nwsLCqL5ramqEuldZaUzOQb7O3LYMAqJtapl6VIYzTskixeWSRUDldT7jyJmI2UjI1v\/rylZdaclOWq7PCPiBAFnRlGJy6656qeY4RKcUTFzIJwRUX+cz5sxZRMRMzj7Ncm6GEXBAgLavKQc0kbXoIVJe1PtMHNGqLePICKQMASbnlEHt3JEdEYtqqa60kEDOYjACJgRou3vrrq\/xav12LP38PXTaYx\/cc+jJjBIjkBYEVF\/nM85y9nLmvHDhQrRr8tbW\/07LLOJOGQFGgBFgBJJGYPv27dE26N8jR47EOgVzKNAgMoqcaUBezpyNM2HChAmgP2F\/2CEs7Bpi+RgBRiBdCMyZMwf0x\/gwOadLG039JnLmPHv2bLRtGzkHI8tZBeuZyTnNE427ZwQYgdAiQNaybj1v2rRJI2om5zSrKxFyVlFpTM5pnmjcPSPACCiBAJ85h0RNTM4hUQSLwQgwAoxACBBgcg6BEnQR+J5ziJTBojACjAAjkEYEmJzTCL6oay\/e2ipua4cMbhaHEWAEGIFQIsDkHEq1OAulutIUhJxFZgQYAUYgpQiovs5nzFUqL1pXWWnsEOZF01yWEWAEshUBldd50hmTM6eMzNZ3l8fNCDACGYwAk7OCylVZaWw5KzjhWGRGgBFIOQIqr\/NsOSsY1o3JOeXvOHfICDACCiLA5MxKSykCTM4phZs7YwQYAUURYHJWUHGqK43C0y1evBgjRoxQIuSocYqw7Ol5YRh3xt0rAirPGRqr6us8O4Qp5hCm+qRT+YVh2b0u7\/6UZ9z9wdFrKyrjrvo6mfVnzsaUkV4nbjrL66nQVJSfZU\/PzGHcGXevCKg8Z2isnDLSq8ZDUJ6Udvnll4OylvDDCDACjAAjkJkI9OnTB48++qiSg8vKbW39q8qYmFtJ7bHQjAAjwAgwArYIqJIKWDSArCVnns+MACPACDACjEBYEWByDqtmWC5GgBFgBBiBrEWAyTlrVc8DZwQYAUaAEQgrAkzOYdUMy8UIMAKMACOQtQgwOWet6nngjAAjwAgwAmFFgMk5rJphuRgBRoARYASyFgEm56xVPQ+cEWAEGAFGIKwIZB05FxYWRnVRU1PjqhcqL1POtSGfCniV36duPTfjVc4w4exFdr1sWOaIrOzGcqRc1eS3TsgwzB8Z7K246+NIN\/4ysuuyeinreeHgClEEsoqcrS+w2wsdxoXX+BK7yZ+uea4yzl5kN5YNgy5kZRfJqpL8ImJO9wdGMtin6z01kq3suiI7znSPKRP6zxpytlt8nH4elq9aksOr\/OmanF7lDNNXuBfZw0ZwKsuezPzWx53Oj4tksU\/Xu+oV97DN+XTiloq+mZxdtq3T+dIbJ4CXBSAVE8euj0TlDAPOXmQP20LlRXaR7tKNfyLyh2Xnwovsxo\/RMHz8e5U93dvv6VzbUt03kzOTs69zzsvLLvPx4atwLo0lKrvd9iDLLo9AItirSs5Wgkvnh5EX3I27FGH4sJCfXWqWzBhyFn2RGieQl0kYNtLwuv2UzqmoMs6Jyu6kn1TpIlHZ9fcm3RaRV\/mt5VUhuLDtWnjBXTRX0ol7qt6tdPWTMeTsBqCXScjk7Iam\/e9VxjkbZQ+D1e\/145OPFBJ\/P601vcz5sOHuHwrhbInJmbe1fZ2ZXl72sH0EJSJ7WCwHL7KHcZH1Kr\/dpE3HDoAX2VW3nMO0Je\/rwhXCxpicmZx9nZaJLlRhIDmvsodBZjfLV5aI0z0Wr9jLWoC+Tm6bxrzILquPVMit+o5FqjBKVz9ZQ86iiSizIMmUSZXywnTO5jTmROQMC86ysodFXqfdh2RJI1Xz2u4DwwvGXsoGMa5k5o0qsie6hgaBdza0mVXkrE8uXbEyWzTpfnFEFoKd\/GGasISbqjjLyG4sY8Q9HduqVoKWwd0qf7rlNhK0jPxhspy9yh5G7GXmvMw4w7QGqS5L1pGz6gpj+RkBRoARYAQyHwEm58zXMY+QEWAEGAFGQDEEmJwVUxiLywgwAowAI5D5CDA5Z76OeYSMACPACDACiiHA5KyYwlhcRoARYAQYgcxHgMk583XMI2QEGAFGgBFQDAEmZ8UUxuIyAowAI8AIZD4CTM6Zr2MeISPACDACjIBiCDA5K6YwFpcRYAQYAUYg8xFgcs58HfMIGQFGgBFgBBRDgMlZMYWxuIwAI8AIMAKZj8D\/A+J45gSfb1sCAAAAAElFTkSuQmCC","height":183,"width":325}}
%---
