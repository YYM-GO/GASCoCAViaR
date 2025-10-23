function [theta,largedelta,smalldelta,nsmall] = threshold2(Delta)
% find optimal theta in "Tail event driven networks of SIFIs" (5)

n = length(Delta);
if n==0
    nn=-1;
else if n>=3
        nmin = ceil(n*0.1);
        nmax = floor(n*0.9);
        nn = nmax-nmin;
    else if n==2
            nmin =1;
            nmax=2;
            nn = nmax-nmin;
        else
            nmin =1;
            nmax=1;
            nn = nmax-nmin;
        end
    end
end


if nn>0
    for i = 1:nn
        ssr(i) = sum((Delta(1:i+nmin)-mean(Delta(1:i+nmin))).^2)+sum((Delta(i+1+nmin:n)-mean(Delta(i+1+nmin:n))).^2);
    end
    nopt = find(ssr == min(ssr));
    theta = (nopt+nmin)/n;
    smalldelta = Delta(1:nopt+nmin);
    largedelta = Delta(nopt+nmin+1:end);
    nsmall = nopt+nmin;
else
    theta=0;
    largedelta=0;
    smalldelta=0;
    nsmall=0;
end

end

