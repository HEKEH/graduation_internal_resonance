function [list_cell, max_size] = cell_to_list_new(sigma2_len, a_cell, idx)
    max_size = 0;
    for k = 1: sigma2_len
        if length(a_cell{idx, k}) > max_size
            max_size = length(a_cell{idx, k});
        end
    end
    list = -100 * ones(max_size, sigma2_len);
    for k = 1: sigma2_len
        temp = a_cell{idx, k};
        list(1: length(temp), k) = temp;
    end
    list_cell = {};
    lc_idx = 1;
    temp_list = [];
    for i = 1: max_size
        last = -100;
        max_interval = 0;
        for j = 1: sigma2_len
            if (list(i,j) ~= -100)% && (last == -100 || max_interval == 0|| abs(list(i,j) - last) < 2 * max_interval))
                if last ~= -100
                    max_interval = max(max_interval, abs(list(i,j) - last));
                end
                temp_list = [temp_list,[j; list(i,j)]];
                last = list(i,j);
            else
                if ~isempty(temp_list)
                    list_cell{lc_idx} = temp_list;
                    lc_idx = lc_idx + 1;
                    temp_list = [];
                    if list(i,j) ~= -100
                        temp_list = [temp_list,[j; list(i,j)]];
                        last = list(i,j);
                    else
                        last = -100;
                    end
                end
            end
            if ~isempty(temp_list)
                    list_cell{lc_idx} = temp_list;
            end
        end
    end
end