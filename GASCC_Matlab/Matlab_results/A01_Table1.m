close
clear
clc


load shenwan       % Shenwan industry index
load hushen        % CSI 300
[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05; %probability level for VaR and CoVaR

[hJBHS,pJBHS] = jbtest(hushen); warning off
[hLBHS,pLBHS] = lbqtest(hushen); warning off
[hLBHS2,pLBHS2] = lbqtest(hushen.^2); warning off
for i = 1:N
    [hJBSW(i),pJBSW(i)] = jbtest(shenwan(:,i));warning off
    [hLBSW(i),pLVSW(i)] = lbqtest(shenwan(:,i));warning off
    [hLBSW2(i),pLVSW2(i)] = lbqtest(shenwan(:,i).^2);warning off
end
sum(hLBSW)
sum(hLBSW2)
sum(hJBSW)


%%                Table 1
% descriptive statistics
% Panel A
A = [mean(hushen),median(hushen),std(hushen),skewness(hushen),kurtosis(hushen)]
% Panel B
for i = 1:N
    a(i,:) = [mean(shenwan(:,i)),median(shenwan(:,i)),std(shenwan(:,i)),skewness(shenwan(:,i)),kurtosis(shenwan(:,i))];
end
B = [mean(a);min(a);quantile(a,0.25);median(a);quantile(a,0.75);max(a)]