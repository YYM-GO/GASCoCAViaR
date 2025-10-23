close
clear
clc
load shenwan
load hushen
load sc1
%---------------- Table 7 ---------------------
daylong = [55 244 244 243 244 243 243 242 242];
cdaylong = [0,cumsum(daylong)];
for k = 1:length(daylong)
    sc1year(k,:) = mean(sc1(cdaylong(k)+1:cdaylong(k+1),:));
end

for k = 1:length(daylong)
    sc1yearsector(k,1) = mean(mean(sc1(cdaylong(k)+1:cdaylong(k+1),1:4)));
    sc1yearsector(k,2) = mean(mean(sc1(cdaylong(k)+1:cdaylong(k+1),5:14)));
    sc1yearsector(k,3) = mean(mean(sc1(cdaylong(k)+1:cdaylong(k+1),15:31)));
end
% results:
sc1year'
mean(sc1year')
sc1yearsector'


