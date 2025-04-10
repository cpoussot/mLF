function ok = check(p_c,p_r)

for ii = 1:length(p_c)
    col_ok  = p_c{ii}.';
    row_ok  = p_r{ii}.';
    col_row = [col_ok; row_ok];
    tmp     = unique(col_row);
    if numel(tmp) == numel(col_row)
        ok = true;
    else
        warning('Interpolation points selection column/row pair may be non unique.')
        ok = false;
        break
    end
end
