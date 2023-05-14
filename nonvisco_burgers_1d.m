% 初期化
clearvars;

% パラメーター
Lx = 4;% 計算領域
dLx = 0.01;% 要素サイズ
dt = 0.005;
nx = Lx / dLx + 1;% 要素数
t_max = 3;
iter_max = t_max / dt;

% 座標の生成
x = 0 : dLx : Lx;

% 配列確保
qf = zeros(nx, 1);% 保存型表示の解
q = zeros(nx, 1);% 非保存型表示の解
flux = zeros(nx - 1, 1);

% 初期条件の設定
for i =  1 : nx
    if abs(x(i) - 1) <= 10^-6
        q(i) = 1;
        qf(i) = 1;
    elseif x(i)  < 1
        q(i) = 1;
        qf(i) = 1;
    else
        q(i) = 0;
        qf(i) = 0;
    end
end

% 時間発展
for iter = 1 : iter_max

    time = iter * dt;
    q_old = q;
    qf_old = qf;

    % 保存型の時間発展
    flux = CalcFlux(flux, qf_old, nx);% フラックスの計算
    for i = 2 : nx - 1
        qf(i) = qf_old(i) - (dt / dLx) * (flux(i) - flux(i - 1));
    end

    % 非保存型の時間発展
    for i = 2 : nx - 1
        q(i) = q_old(i)...
            - 0.5 * (q_old(i) + abs(q_old(i))) * (dt / dLx) * (q_old(i) - q_old(i - 1))...
            - 0.5 * (q_old(i) - abs(q_old(i))) * (dt / dLx) * (q_old(i + 1) - q_old(i));
    end

    % コマンドウィンドウへの出力
    txt = ['iter = ', num2str(iter), ' / ', num2str(iter_max)];
    disp(txt);

    % 動画保存
    if iter == 1
        plotconfig(x, q, qf, time)
        filename = ['nonvisco_burgers_1d','.mp4'];
        v = VideoWriter(filename, 'MPEG-4');
        v.FrameRate = 40;
        open(v);
    else
        plotconfig(x, q, qf, time)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

end

% 動画ファイルを閉じる
close(v);

function [flux] = CalcFlux(flux, q, nx)
% 保存型表示のフラックスを求める。
% flux(i)はf_{i+1/2}に相当

for i = 1 : nx - 1

    ur = q(i + 1);
    ul = q(i);
    fr = 0.5 * ur^2;
    fl = 0.5 * ul^2;
    c = 0.5 * (ur + ul);
    flux(i) = 0.5 * (fr + fl -sign(c) * (fr - fl));

end

end

function [] = plotconfig(x, qf, q, t)

plot(x, qf, x, q)

title(['time = ', num2str(t, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
axis equal;
axis tight;
axis on;
fig=gcf;
fig.Color='white';
ylim([-0.2 1.4]);
xlabel('x')
ylabel('q')

% 凡例
legend({'non-conservative', 'conservative'},'Location','southeast','FontSize', 10)

% 新しいプロットの時、軸設定を保持したまま前のグラフィックスオブジェクトを消去
set(gca,'nextplot','replacechildren');

end

