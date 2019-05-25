%本文件的作用是设置索梁结构的物理参数。
%参数默认单位: 质量: Kg, 尺寸: m, 时间: s, 力: N;
personal_setting;
g = 9.81;                               % 重力加速度
theta = 30 / 180 * pi;            % 索、梁初始夹角
Ec = 200 * 1e9;                     %索弹性模量
Ac = 3.14 * 10 ^ (-2);            %索面积
rhoc = 7850;                        %索密度
mc = rhoc * Ac;                     %索的单位长度质量
Hc = 2e7;                               %索的初始张力
lc = 100;                               %索的长度 
H = Hc / (mc * g * lc * cos(theta));  %索的无量纲初始索力

%----------------------计算斜拉索的Irvine常数---------------------------
zc = @(xc)1/2 * (mc * g * cos(theta) / Hc .* xc .* (lc - xc));
zc_prime = @(xc)(-mc * g / Hc * cos(theta) .* (xc - lc / 2));
s_prime_3 = @(xc)(sqrt(1 + zc_prime(xc) .^ 2)) .^ 3;
Le = integral(s_prime_3, 0, lc);
lambda = sqrt((mc * g * cos(theta) * lc / Hc)^2 * (Ec * Ac * lc) / (Hc * Le)); %Irvine常数
%------------------------------------------------------------------------

Eb = 3.45e10;                     %梁的弹性模量
lb = lc * cos(theta);            %梁的长度
mb = 4 * 10^4;                 %梁的单位长度质量
Ib = K * Ec * Ac^2 / Eb;    %梁的惯性矩

%------------------------其余无量纲化参数计算-------------------
m = mc / mb;
P = mc / mb / cos(theta);
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 = m / (beta1 * cos(theta));
k2 = m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 = m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 = -m / (beta1 * cos(theta));
%-----------------------------------------------------------------
