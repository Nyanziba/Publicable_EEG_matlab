%[text] # subplotの使い方
% 関数 subplot -> タイル状に配置された Axes の作成
figure
for i=1:6
    subplot(2,3,i)
end

% でも，これだと決まった並び方しかできないのでカッコ悪い。
% subplot('Position',...を使う。

% subplot('Position',positionVector) は、positionVector で指定された位置に
% 新しい座標軸を作成します。positionVector は、エントリが 0.0 ～ 1.0 の範囲の
% 正規化された値になる [left,bottom,width,height] の形式の 4 要素ベクトルです。
% 位置ベクトルで以前の Axes をオーバーラップする Axes を指定した場合、既存の 
% Axes が新しい Axes によって置き換えられます。
%%
%[text] ## plot ERP data(individual)
fs = 1000; % サンプリング周波数
N=size(data,2); % データ長

% ----------trigger
marker=zeros(2,2);
marker(:,1)=0.0;
marker(1,2)=-30;
marker(2,2)=30;
% ----------trigger

% x=linspace(0,N/fs,N); 
x=linspace(-0.1,(N-100)/fs,N);

load('EEG_locs_new.mat'); % channel位置を取得する。

sz = get(0, 'ScreenSize'); % スクリーンサイズを取得する。


% data->(ch,time,epoch,dev,prt); (31,800,1000,3,X);

for prt=1:size(data,5)
    figure('Position', [sz(1) sz(2) sz(3) sz(4)]) % 画面いっぱいのfigureを作成する
    for chl=1:size(sloc,1)
        subplot('position',[sloc(chl,1) sloc(chl,2) 0.07 0.07])
        hold on
        % ---------------------------- y1-y6をMMNanalysis004.mlxからコピペ
        %findなんてしないで、論理演算子で
        y1=mean(data(chl, :, trg(:,1,prt) == 1, 1, prt),3 ,"omitnan"); % 1010Hz % standard
        y2=mean(data(chl, :, trg(:,1,prt) == 2, 1, prt),3 ,"omitnan"); % 1010Hz % deviant
        y3=mean(data(chl, :, trg(:,2,prt) == 1, 2, prt),3 ,"omitnan"); % 1030Hz % standard
        y4=mean(data(chl, :, trg(:,2,prt) == 2, 2, prt),3 ,"omitnan"); % 1030Hz % deviant
        y5=mean(data(chl, :, trg(:,3,prt) == 1, 3, prt),3 ,"omitnan"); % 1050Hz % standard
        y6=mean(data(chl, :, trg(:,3,prt) == 2, 3, prt),3 ,"omitnan"); % 1050Hz % deviant
        % ----------------------------
        plot(x,y1,'Color', [255  75   0]/255,'LineWidth',2,'LineStyle','--'); % 1000Hz % standard
        plot(x,y2,'Color', [255  75   0]/255,'LineWidth',2); % 1010Hz % deviant
        plot(x,y3,'Color', [  0  90 255]/255,'LineWidth',2,'LineStyle',':'); % 1000Hz % standard
        plot(x,y4,'Color', [  0  90 255]/255,'LineWidth',2); % 1030Hz % deviant
        plot(x,y5,'Color', [  3 175 122]/255,'LineWidth',2,'LineStyle','-.'); % 1000Hz % standard
        plot(x,y6,'Color', [  3 175 122]/255,'LineWidth',2); % 1050Hz % deviant
        axis('ij')
        ymx=25;
        xlim([-0.1 0.7])
        ylim([-ymx ymx])
        title([ ch_name{chl} '(' num2str(chl) ')' ]) % titleを付けます。
        plot(marker(:,1),marker(:,2),'k:','LineWidth',2); % trigger
        set(gca,'box','on');
        set(gca,'linewidth',1)
    end
    legend('1000Hz', '1010HzDev', '1000Hz', '1030HzDev', '1000Hz', '1050HzDev')
    makeSubplotsClickable();
    saveas(gcf,[fn0{prt}],'fig'); % 拡張子.figでfigureを保存
    % saveas(gcf,[fn0{par}],'emf'); % 拡張子.emfでfigureを保存
    % close(gcf);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":25.2}
%---
