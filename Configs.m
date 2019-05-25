%���ļ������������������ṹ�����������
%����Ĭ�ϵ�λ: ����: Kg, �ߴ�: m, ʱ��: s, ��: N;
personal_setting;
g = 9.81;                               % �������ٶ�
theta = 30 / 180 * pi;            % ��������ʼ�н�
Ec = 200 * 1e9;                     %������ģ��
Ac = 3.14 * 10 ^ (-2);            %�����
rhoc = 7850;                        %���ܶ�
mc = rhoc * Ac;                     %���ĵ�λ��������
Hc = 2e7;                               %���ĳ�ʼ����
lc = 100;                               %���ĳ��� 
H = Hc / (mc * g * lc * cos(theta));  %���������ٳ�ʼ����

%----------------------����б������Irvine����---------------------------
zc = @(xc)1/2 * (mc * g * cos(theta) / Hc .* xc .* (lc - xc));
zc_prime = @(xc)(-mc * g / Hc * cos(theta) .* (xc - lc / 2));
s_prime_3 = @(xc)(sqrt(1 + zc_prime(xc) .^ 2)) .^ 3;
Le = integral(s_prime_3, 0, lc);
lambda = sqrt((mc * g * cos(theta) * lc / Hc)^2 * (Ec * Ac * lc) / (Hc * Le)); %Irvine����
%------------------------------------------------------------------------

Eb = 3.45e10;                     %���ĵ���ģ��
lb = lc * cos(theta);            %���ĳ���
mb = 4 * 10^4;                 %���ĵ�λ��������
Ib = K * Ec * Ac^2 / Eb;    %���Ĺ��Ծ�

%------------------------���������ٻ���������-------------------
m = mc / mb;
P = mc / mb / cos(theta);
beta1 = Eb * Ib * lc^2 * mc / (Hc * lb^4 * mb);
k1 = m / (beta1 * cos(theta));
k2 = m * (2 * H * tan(theta) - 1) * lambda^2 / (2 * beta1 * cos(theta));
k3 = m * (- 2 * H * tan(theta) + 1)^2 * lambda^2 * cos(theta) / (4 * beta1);
k4 = -m / (beta1 * cos(theta));
%-----------------------------------------------------------------
