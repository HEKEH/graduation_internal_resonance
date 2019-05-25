function list_cell_i = merge_items(list_cell_i, idx1, idx2)
    list_cell_i{idx1} = [list_cell_i{idx1},list_cell_i{idx2}];
    list_cell_i{idx2} = [];
end