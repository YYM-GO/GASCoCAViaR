function out = findposition2(position)
% find the position (lower and right) in the triangle using the vec index.
N=31;
out = zeros(length(position),2);
collength = [N-1:-1:1];
collengthall = [0,cumsum(collength)];
for k = 1:length(position)
   j =  max(find((collengthall-position(k))<0));       % in which column
   i = position(k)-collengthall(j)+j;                  % in which row
   out(k,:) = [i,j];
end


