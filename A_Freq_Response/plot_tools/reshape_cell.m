function list_cell = reshape_cell(list_cell)
    for idx = 1: length(list_cell)
        temp_cell = list_cell{idx};
        temp_idx = [];
        for j = 1: length(temp_cell)
            if ~(isempty(temp_cell{j}) || length(temp_cell{j}) <= 4)
                temp_idx = [temp_idx, j];
            end
        end
        list_cell{idx} = list_cell{idx}(temp_idx);
    end
end

