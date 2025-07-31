% Matrix size estimation in Mo
% with double precision
% input is the dimensions of the matrix

function s = file_size(mat,type)

if nargin < 2
    type = 'real';
end
%
s = prod(mat)*8/(2^10^2);
%
if strcmp(type,'complex')
    s = 2*s; % double complex
end

if s/(2^10)^2 > 1
    fprintf(['Space on disk (' type '): %3.2f To\n'],s/(2^10)^2)
elseif s/(2^10) > 1
    fprintf(['Space on disk (' type '): %3.2f Go\n'],s/(2^10))
elseif s >= 1
    fprintf(['Space on disk (' type '): %3.2f Mo\n'],s)
elseif s < 1
    fprintf(['Space on disk (' type '): %3.2f Ko\n'],s*(2^10))
end