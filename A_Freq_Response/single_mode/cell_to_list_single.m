function [list, max_size, start_idx, end_idx, left] = cell_to_list_single(sigma_len, some_cell, left)
    max_size = 0;
    for k = 1: sigma_len
        if length(some_cell{k}) > max_size
%             idx
%             k
%             some_cell{k}
            max_size = length(some_cell{k});
        end
    end
    list = zeros(max_size, sigma_len);
    start_idx = 0;
    end_idx = 0;
    %max_size
    if max_size == 1
        for k = 1: sigma_len
            if isempty(some_cell{k})
                list(k) = 0;
            else
                list(k) = some_cell{k};
            end
        end
        left = 0;
    elseif max_size == 3
        for k = 1: sigma_len
            if start_idx == 0 
                if length(some_cell{k}) > 1
                    start_idx = k;
                end
            else 
                if length(some_cell{k}) == 1
                    end_idx = k - 1;
                    break;
                end
            end
        end
        if (nargin < 3)
            if abs(max(some_cell{start_idx}) - some_cell{start_idx - 1}) > abs(min(some_cell{start_idx}) - some_cell{start_idx - 1})
                left = 1;
            else
                left = 0;          
            end
        end
        for k = 1: end_idx
            tmpk = some_cell{k};
            if left
                list(1, k) = tmpk(1);
            else
                list(1, k) = tmpk(end);
            end
        end
        if left == 1
            list(2, start_idx) = some_cell{start_idx}(end);
        else
            list(2, start_idx) = some_cell{start_idx}(1);
        end
        for k = start_idx + 1: end_idx
            tmpk = some_cell{k};
            list(2, k) = tmpk(2);
        end
        if left == 1
            list(2, end_idx) = some_cell{end_idx}(1);
        else
            list(2, end_idx) = some_cell{end_idx}(end);
        end
        for k = start_idx: sigma_len
            tmpk = some_cell{k};
            if left == 1
                list(3, k) = tmpk(end);
            else
                list(3, k) = tmpk(1);
            end
        end
    else
        left = 0;
        '´íÎó,²»µÈÓÚ3';
    end
end