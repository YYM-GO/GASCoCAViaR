close
clear
clc
load shenwan
load hushen
load S


%%------------ Figure 4 ----------------------------
daylong = [55 244 244 243 244 243 243 242 242];
cdaylong = [0,cumsum(daylong)];
y = [{' '},{'2016'},{'2017'},{'2018'},{'2019'},{'2020'},{'2021'},{'2022'},{'2023'},{'2024'}];
dastr=[]
for i = 1:length(y)
    dastr{i} = y{i};
end
datestrvec = reshape(dastr',[],1);
dastr1 = datestrvec(end-length(cdaylong)+1:end);
dastr2 = cell(length(cdaylong),1);
for i = 1:length(cdaylong)
    ymlength2(i) = cdaylong(i);
    dastr2{i} = dastr1{i};
end

plot(S,'LineWidth',0.9,'Color','#0072BD'),title('Systemic Risk Score');
set(gca,'XTick',cdaylong);
set(gca,'XTickLabel',dastr1)
xlim([0 cdaylong(end)])
set(gca,'XGrid','off')


