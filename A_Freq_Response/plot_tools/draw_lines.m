function draw_lines(sigma2_range, a_cell)
    lgd = {'\fontname{Times new Roman}{\ita}_m', '\fontname{Times New Roman}{\ita}_n'};
    colors = {'r', [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],...
        [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
    
    sigma2_len = length(sigma2_range);
    list_cell = cell(2,1);
    for idx = 1: 2
        [list_cell{idx}, max_size] = cell_to_list_new(sigma2_len, a_cell, idx);
    end 

    list_cell_size = length(list_cell{1}) + length(list_cell{2});
    p_idx = 1;
    test_plot;
    list_cell{1} = divide(list_cell{1}, 1);
    list_cell{1} = divide(list_cell{1}, 1);
    list_cell = reshape_cell(list_cell);
    list_cell{1} = add_element(list_cell{1}, 2, 6, '开头', '开头');
    list_cell{1} = add_element(list_cell{1}, 3, 5, '结束', '结束');
    
    list_cell{2} = divide(list_cell{2}, 1);
    list_cell{2} = divide(list_cell{2}, 1);
    list_cell = reshape_cell(list_cell);
    list_cell{2} = add_element(list_cell{2}, 2, 6, '开头', '开头');
    test_plot;

    lines = cell(2, 1);
    for i = 1: 2
        temp_cell = list_cell{idx};
        for j = 1: length(temp_cell)
            lines{i}{j} = '-';
        end
    end
    lines{1}{2} = '--';
    lines{1}{5} = '--';
    lines{2}{2} = '--';
    lines{2}{5} = '--';

    for idx = 1: 2
        temp_cell = list_cell{idx};
        figure(p_idx)
        p_idx = p_idx + 1;
        linewidth = 2;
        for j = 1: length(temp_cell)
            xx = sigma2_range(temp_cell{j}(1,:));
            yy = temp_cell{j}(2,:);
            PL = plot(xx, yy, lines{idx}{j} , 'Color', colors{idx},'linewidth', linewidth);
            hold on
        end
        hold off
        set(gca,'FontSize',20)
        set(gcf,'Position',[50,50,1000,750], 'color','w') 
        set(gca,'Position',[0.1 0.12 0.85 0.8]);
        xlabel('\fontname{Times New Roman}{\it\sigma}_2', 'FontSize', 20)
        ylabel(lgd{idx}, 'FontSize', 20)
        xlim([sigma2_range(1), sigma2_range(end)]);
        ylim([0, Inf]);
    end
end