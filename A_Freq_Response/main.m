clear;
addpath('func_tools');
addpath('plot_tools');
choice = 'n';
run('get_p.m');
sigma2_range = -6e-2: 5e-5: 6e-2;
sigma2_len = length(sigma2_range);
a_cell = cell(2, sigma2_len);%第一行m，第二行n
gamma1_angle = cell(1, sigma2_len);

if choice == 'm'
    get_a_n = @(a_m, sigma2) [-((a_m.^2.*pn(2)+4.*sigma1.*omega_n+4.*sigma2.*omega_n-2.*sqrt(a_m.^2.*(m.*J_n+pn(4).*sigma1).^2-4.*((-1+pn(1)).*mu_b-pn(1).*mu_c).^2.*omega_n.^2))./pn(3));...
              -((a_m.^2.*pn(2)+4.*sigma1.*omega_n+4.*sigma2.*omega_n+2.*sqrt(a_m.^2.*(m.*J_n+pn(4).*sigma1).^2-4.*((-1+pn(1)).*mu_b-pn(1).*mu_c).^2.*omega_n.^2))./pn(3))];
    get_ret = @(a_m, a_n, sigma2) -a_m.^2.*(pm(5)+pm(6)).^2.*(m.*J_n+pn(4).*sigma1).^2+1./4.*(a_m.^2.*(((-(((-1)+pm(1)))).*mu_b+pm(1).*mu_c)).*((m.*J_n+pn(4).*sigma1)).*omega_m+a_n.^2.*(((-(((-1)+pn(1)))).*mu_b+pn(1).*mu_c)).*((m.*J_m+pm(4).*sigma1)).*omega_n).^2+1./64.*(a_n.^2.*((a_m.^2.*pn(2)+a_n.^2.*pn(3))).*((m.*J_m+pm(4).*sigma1))-a_m.^2.*((a_m.^2.*pm(2)+a_n.^2.*pm(3))).*((m.*J_n+pn(4).*sigma1))-8.*a_m.^2.*(m.*J_n+pn(4).*sigma1).*sigma2.*omega_m+4.*a_n.^2.*(m.*J_m+pm(4).*sigma1).*(sigma1+sigma2).*omega_n).^2;
    get_cos_gamma1 = @(a_m, a_n, sigma2) -((1./2.*(a_m.^2.*pn(2)+a_n.^2.*pn(3))+2.*omega_n(sigma1+sigma2))./a_m(m.*J_n+pn(4).*sigma1));
    get_sin_gamma1 = @(a_m, a_n, sigma2) (2*((1-pn(1)).*mu_b+pn(1).*mu_c).*omega_n)./(a_m.*(m.*J_n+pn(4).*sigma1));
else
    get_a_n = @(a_m, sigma2)sqrt([(a_m.^2.*pm(3).*((a_m.^2.*pm(2)-8.*((sigma1-2.*sigma2)).*omega_m))+2.*sqrt(a_m.^6.*pm(2).^2.*(m.*J_m+pm(4).*sigma1).^2+16.*a_m.^2.*(m.*J_m+pm(4).*sigma1).^2.*(((-1+pm(1)).*mu_b-pm(1).*mu_c).^2+4.*(sigma1-2.*sigma2).^2).*omega_m.^2+4.*a_m.^4.*omega_m.*(-4.*pm(2).*(m.*J_m+pm(4).*sigma1).^2.*(sigma1-2.*sigma2)-pm(3).^2.*((((-1)+pm(1))).*mu_b-pm(1).*mu_c).^2.*omega_m)))./((2.*m.*J_m-a_m.*pm(3)+2.*pm(4).*sigma1).*(2.*m.*J_m+a_m.*pm(3)+2.*pm(4).*sigma1));...
              (a_m.^2.*pm(3).*((a_m.^2.*pm(2)-8.*((sigma1-2.*sigma2)).*omega_m))-2.*sqrt(a_m.^6.*pm(2).^2.*(m.*J_m+pm(4).*sigma1).^2+16.*a_m.^2.*(m.*J_m+pm(4).*sigma1).^2.*(((-1+pm(1)).*mu_b-pm(1).*mu_c).^2+4.*(sigma1-2.*sigma2).^2).*omega_m.^2+4.*a_m.^4.*omega_m.*(-4.*pm(2).*(m.*J_m+pm(4).*sigma1).^2.*(sigma1-2.*sigma2)-pm(3).^2.*((((-1)+pm(1))).*mu_b-pm(1).*mu_c).^2.*omega_m)))./((2.*m.*J_m-a_m.*pm(3)+2.*pm(4).*sigma1).*(2.*m.*J_m+a_m.*pm(3)+2.*pm(4).*sigma1))]);
    get_ret = @(a_m, a_n, sigma2) -a_n.^2.*(pn(5)+pn(6)).^2.*(m.*J_m+pm(4).*sigma1).^2+1./4.*(a_m.^2.*(((((-1)+pm(1))).*mu_b-pm(1).*mu_c)).*((m.*J_n+pn(4).*sigma1)).*omega_m-a_n.^2.*(((((-1)+pn(1))).*mu_b-pn(1).*mu_c)).*((m.*J_m+pm(4).*sigma1)).*omega_n).^2+1./64.*(a_m.^2.*((m.*J_n+pn(4).*sigma1)).*((a_m.^2.*pm(2)+a_n.^2.*pm(3)-8.*((sigma1-2.*sigma2)).*omega_m))-a_n.^2.*((m.*J_m+pm(4).*sigma1)).*((a_m.^2.*pn(2)+a_n.^2.*pn(3)+8.*sigma2.*omega_n))).^2;
    get_cos_gamma1 = @(a_m, a_n, sigma2)-((a_m.*(1./2.*(a_m.^2.*pm(2)+a_n.^2.*pm(3))+4.*(-sigma1+2.*sigma2).*omega_m))./(a_n.^2.*((m.*J_m+pm(4).*sigma1))));
    get_sin_gamma1 = @(a_m, a_n, sigma2)-((2*((1-pm(1)).*mu_b+pm(1).*mu_c).*a_m.*omega_m)./(a_n.^2.*((m.*J_m+pm(4).*sigma1))));
end

valid_idx = [];
for sidx = 1: sigma2_len
    sigma2 = sigma2_range(sidx);
    [a_m,a_n] = get_a(get_a_n, get_ret, sigma2, @get_ret_2);
    if ~isempty(a_m)
        valid_idx = [valid_idx, sidx];
        a_cell{1, sidx} = a_m;
        a_cell{2, sidx} = a_n;
        
        tmp = get_cos_gamma1(a_m, a_n, sigma2);
        tmp(tmp < -1) = -1;
        tmp(tmp > 1) = 1;
        angle = acos(tmp);
        nn_idx = get_sin_gamma1(a_m, a_n, sigma2) < 0;
        angle(nn_idx) =  - angle(nn_idx);
        gamma1_angle{1, sidx} = angle;
    end
end
if isempty(valid_idx)
    disp('无有效数据');
else
    sigma2_range = sigma2_range(valid_idx);
    a_cell = a_cell(:, valid_idx);
    sigma2_len = length(sigma2_range);
    colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
        [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
    
    draw_lines(sigma2_range, a_cell);
    
    am_an_bar;
    draw_bar_lines(sigma2_range, a_bar_cell);
    
end