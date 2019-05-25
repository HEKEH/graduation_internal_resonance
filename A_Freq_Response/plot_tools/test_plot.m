for idx = 1: 2
    line_idx = 1;
    line_lgd = {};
    PLine = [];
    figure(p_idx)
    p_idx = p_idx + 1;
    linewidth = 2;
    temp_cell = list_cell{idx};
    temp_idx = [];
    for j = 1: length(temp_cell)
        if ~(isempty(temp_cell{j}) || length(temp_cell{j}) <= 4)
            temp_idx = [temp_idx, j];
        end
    end
    list_cell{idx} = list_cell{idx}(temp_idx);
    temp_cell = list_cell{idx};
    for j = 1: length(temp_cell)
        temp_color = rand(1,3);
        xx = sigma2_range(temp_cell{j}(1,:));
        yy = temp_cell{j}(2,:);
        temp_len = length(xx);
        PLine(line_idx) = plot(xx, yy, 'Color',temp_color);
        line_lgd{line_idx} = num2str(j);
        text(xx(ceil(temp_len / 2)), yy(ceil(temp_len / 2)), num2str(j), 'Color',temp_color);
        line_idx = line_idx + 1;
        hold on
    end
    hold off
    legend(PLine, line_lgd, 'location', 'northwestoutside')
    set(gca,'FontSize',14)
    set (gcf,'Position',[50,50,1000,750], 'color','w') 
    xlabel('\it\sigma_2', 'FontSize', 15)
    ylabel(['\it', lgd{idx}], 'FontSize', 15)
    ylim([0, Inf]);
end