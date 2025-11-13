function p_rnd = rand(n,lam,COMPLEX)

p_rnd = [];
%n = numel(lam);
if ischar(COMPLEX)
    if strcmp(COMPLEX,'R')
        COMPLEX = false;
    end
    if strcmp(COMPLEX,'C')
        COMPLEX = true;
    end
end

COMPLEX_INIT = COMPLEX;
if iscell(lam)
    for ii = 1:n
        COMPLEX = 0;
        REAL    = 1;
        if ii == 1
            COMPLEX = COMPLEX_INIT;
            REAL    = not(COMPLEX_INIT);
        end
        a           = min(lam{ii});
        b           = max(lam{ii});
        rnd_r       = (b-a).*rand(1) + a;
        rnd_i       = (b-a).*rand(1) + a;
        p_rnd(ii)   = rnd_r*REAL+1i*rnd_i*COMPLEX;
    end
else
    a       = min(lam);
    b       = max(lam);
    rnd_r   = (b-a).*rand(1,n) + a;
    rnd_i   = (b-a).*rand(1,n) + a;
    p_rnd   = rnd_r+1i*rnd_i*COMPLEX;
end