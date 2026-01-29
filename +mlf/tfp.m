function Htf = tfp(p_c,w,c,param)

n = length(p_c);
for ii = 1:n
    dim(ii) = length(p_c{ii});
end
comb    = mlf.combinations_dim(dim);
c       = double(c);
w       = double(w);
N       = length(c);
Nstep   = prod(dim(2:end));

for ii = 1:N
    tmp = num2cell(comb(ii,:));
    for jj = 1:n
        ip(ii,jj) = p_c{jj}(tmp{jj});
    end
end
for ii = 1:length(c)
    p_i     = ip(ii,2:end);
    d(ii,1) = 1/prod(param-p_i);
end
c   = c.*d;
idx = find(imag(ip(:,1))>0);
c   = c(idx);
w   = w(idx);
ip  = ip(idx,:);
N2  = length(p_c{1})/2;

%syms s
for ii = 1:length(c)
    c_i = c(ii);
    w_i = w(ii);
    s_i = ip(ii,1);
    x   = real(c_i);
    y   = imag(c_i);
    a   = real(s_i);
    b   = imag(s_i);
    % DEN
    d_n(ii,:) = [2*x -2*a*x-2*b*y];
    d_d(ii,:) = [1 -2*a a^2+b^2];
    %dd(ii,1)  = c_i/(s-s_i) + conj(c_i)/(s-conj(s_i));
    % NUM
    x   = real(w_i*c_i);
    y   = imag(w_i*c_i);
    a   = real(s_i);
    b   = imag(s_i);
    n_n(ii,:) = [2*x -2*a*x-2*b*y];
    n_d(ii,:) = [1 -2*a a^2+b^2];
    %nn(ii,1)  = w_i*c_i/(s-s_i) + conj(w_i*c_i)/(s-conj(s_i));
end
%
NUM_ = 0;
DEN_ = 0;
offset = 0;
for ii = 1:N2
    range = offset+(1:Nstep);
    offset = range(end);
    NUM_ = NUM_ + tf(sum(n_n(range,:),1),n_d(range(1),:));
    DEN_ = DEN_ + tf(sum(d_n(range,:),1),d_d(range(1),:));
end
Hss = minreal(ss(NUM_/DEN_));
Htf = tf(Hss);
