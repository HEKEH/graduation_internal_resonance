function list_cell_i = divide(list_cell_i, idx1)
    list = list_cell_i{idx1};
    [~, qd_idx] = max(abs(list(2, 2: end)-list(2, 1: end-1)));
    list_cell_i{idx1} = list(:, 1: qd_idx);
    list_cell_i{length(list_cell_i) + 1} = list(:, qd_idx + 1: end);
end