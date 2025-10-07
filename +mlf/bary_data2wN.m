function [c_new,curr_new] = bary_data2wN(c,curr)

c_new       = c;
curr_new    = curr;
%
field_name1 = fieldnames(c);
field_name2 = fieldnames(curr{1});
if strcmp(field_name2{1},'nflop')
    n       = str2double(field_name2{2}(2:end));
else
    n       = str2double(field_name2{1}(2:end));
end
for ii = 1:numel(field_name1)
    eval(['c' num2str(ii) '= c.' field_name1{ii} ';']);
end
for jj = 1:length(curr)
    eval(['c_new.w' num2str(n) ' = [c_new.w' num2str(n) ' curr{jj}.w' num2str(n) '];']);
    if isfield(curr{jj},'i2')
        [c1,c2]     = mlf.bary_data2w12(curr{jj},c1,c2);
        c_new.w1    = c1;
        c_new.w2    = c2;
    else
        eval(['curr_new = curr{' num2str(jj) '}.i' num2str(n-1) ';']);
        [c_new,curr_new] = mlf.bary_data2wN(c_new,curr_new);
    end
end