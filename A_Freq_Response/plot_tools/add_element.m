function list_cell_i = add_element(list_cell_i, idx1, idx2, l1, l2)
    if l2 == '开头'
        elem = list_cell_i{idx2}(:, 1);
    else
        elem = list_cell_i{idx2}(:, end);
    end
    if l1 == '开头'
        list_cell_i{idx1} = [elem, list_cell_i{idx1}];
    else
        list_cell_i{idx1} = [list_cell_i{idx1}, elem];
    end
end