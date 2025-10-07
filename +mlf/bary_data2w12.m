function [c1,c2] = bary_data2w12(data,c1,c2)

tmp = data.i2;
for kk = 1:numel(tmp)
    c2 = [c2 tmp{kk}.w2];
    c1 = [c1 tmp{kk}.w1];
end
