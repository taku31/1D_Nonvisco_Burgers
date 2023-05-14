% ������
clearvars;

% �p�����[�^�[
Lx = 4;% �v�Z�̈�
dLx = 0.01;% �v�f�T�C�Y
dt = 0.005;
nx = Lx / dLx + 1;% �v�f��
t_max = 3;
iter_max = t_max / dt;

% ���W�̐���
x = 0 : dLx : Lx;

% �z��m��
qf = zeros(nx, 1);% �ۑ��^�\���̉�
q = zeros(nx, 1);% ��ۑ��^�\���̉�
flux = zeros(nx - 1, 1);

% ���������̐ݒ�
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

% ���Ԕ��W
for iter = 1 : iter_max

    time = iter * dt;
    q_old = q;
    qf_old = qf;

    % �ۑ��^�̎��Ԕ��W
    flux = CalcFlux(flux, qf_old, nx);% �t���b�N�X�̌v�Z
    for i = 2 : nx - 1
        qf(i) = qf_old(i) - (dt / dLx) * (flux(i) - flux(i - 1));
    end

    % ��ۑ��^�̎��Ԕ��W
    for i = 2 : nx - 1
        q(i) = q_old(i)...
            - 0.5 * (q_old(i) + abs(q_old(i))) * (dt / dLx) * (q_old(i) - q_old(i - 1))...
            - 0.5 * (q_old(i) - abs(q_old(i))) * (dt / dLx) * (q_old(i + 1) - q_old(i));
    end

    % �R�}���h�E�B���h�E�ւ̏o��
    txt = ['iter = ', num2str(iter), ' / ', num2str(iter_max)];
    disp(txt);

    % ����ۑ�
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

% ����t�@�C�������
close(v);

function [flux] = CalcFlux(flux, q, nx)
% �ۑ��^�\���̃t���b�N�X�����߂�B
% flux(i)��f_{i+1/2}�ɑ���

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

% �}��
legend({'non-conservative', 'conservative'},'Location','southeast','FontSize', 10)

% �V�����v���b�g�̎��A���ݒ��ێ������܂ܑO�̃O���t�B�b�N�X�I�u�W�F�N�g������
set(gca,'nextplot','replacechildren');

end

