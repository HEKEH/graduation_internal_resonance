clear;
addpath('single_mode')
choice = 'm';
run('./get_p.m');
sigma2_range = -2e-3: 1e-6: 2e-3;
sigma2_len = length(sigma2_range);
a_cell = cell(1, sigma2_len);%第一行m，第二行n
gamma1_angle = cell(1, sigma2_len);

temp1 = 1/ 64 * pm(2) ^ 2;
temp4 = (pm(5) + pm(6))^2;

for sidx = 1: sigma2_len
    sigma2 = sigma2_range(sidx);
    a = get_a_single(pm, temp1, temp4, sigma2, mu_c, mu_b, omega_m);
    a_cell{sidx} = a;
    tmp = (-1/8 * pm(2) * a.^3 - omega_m * sigma2 * a) / (pm(5) + pm(6));
    tmpm(tmp < -1) = -1;
    tmpm(tmp > 1) = 1;
    angle = acos(tmp);
    n_idx = (1/ 2 * omega_m * (pm(1) * mu_c + (1 - pm(1)) * mu_b) * a) / (pm(5) + pm(6)) < 0;
    angle(n_idx) = -angle(n_idx);
    gamma1_angle{sidx} = angle;
end

colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
figure(1)
hold off
linewidth = 2;

[list, max_size, start_idx, end_idx, left] = cell_to_list_single(sigma2_len, a_cell); %

lgd = {'\fontname{Times new Roman}{\ita}_m','\fontname{Times new Roman}{\ita}_n'};
PL = zeros(1, 2);
if max_size == 1
    [~, max_idxs] = max(list);
    PL(1) = plot(sigma2_range, list, 'color', colors{1} ,'LineWidth', linewidth);
else
    [~, max_idxs] = max(max(list));
    PL(1) = plot(sigma2_range(1: end_idx), list(1, 1: end_idx), 'color', colors{1} ,'LineWidth', linewidth);
    hold on
    plot(sigma2_range(start_idx: end_idx), list(2, start_idx: end_idx), ':', 'color', colors{1}, 'LineWidth', linewidth);
    hold on
    plot(sigma2_range(start_idx: sigma2_len), list(3, start_idx: sigma2_len), 'color', colors{1} ,'LineWidth', linewidth);
end
hold on
PL(2) = plot(sigma2_range, zeros(1, length(sigma2_range)), 'color', colors{2} ,'LineWidth', linewidth);
legend(PL, lgd, 'FontSize', 15);

hold off
set(gca,'FontSize',15)
set (gcf,'Position',[50,50,1000,750], 'color','w') 
set(gca,'Position',[0.075 0.1 0.85 0.85]);
xlabel('\fontname{Times new Roman}{\it\sigma}_2', 'FontSize', 15)
ylabel(['\fontname{Times new Roman}\ita', ''], 'FontSize', 15)
xlim([sigma2_range(1), sigma2_range(end)]);
ylim([0, Inf]);