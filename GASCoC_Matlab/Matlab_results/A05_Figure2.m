close
clear
clc
%%%%         plot
load shenwan       % Shenwan industry index
load hushen        % CSI 300
load VaRPOTi2s

[T,N] = size(shenwan);
THETA1 = 0.05; THETA2 = 0.05;

load CoESfullPOTASi2s
load CoVaRfullPOTASi2s
load CoVaRfullPOTASi2sc
load CoESfullPOTASi2sc

load CoESfullPOTASs2i
load CoVaRfullPOTASs2i
load CoVaRfullPOTASs2ic
load CoESfullPOTASs2ic

load ScorePOTASi2s
load gammafullPOTASi2s




y = [{'2010'},{'2011'},{'2012'},{'2013'},{'2014'},{'2015'},{'2016'},{'2017'},{'2018'},{'2019'},{'2020'},{'2021'},{'2022'},{'2023'},{'2024'}];
dastr=[]
for i = 1:length(y)
    dastr{i} = y{i};
end
datestrvec = reshape(dastr',[],1);
daylong = [241 244 243 238 245 244 244 244 243 244 243 243 242 242];

ymlength1 = [0,cumsum(daylong)];
dastr1 = datestrvec(end-length(ymlength1)+1:end);

dastr2 = cell(length(ymlength1),1);
for i = 1:length(ymlength1)
    ymlength2(i) = ymlength1(i);
    dastr2{i} = dastr1{i};
end
%%%%%%%%%%%%%%%%%%%%%% figure in the main text %%%%%%%%%%%%%%%%%%%%%%%%
[ha, pos] = tight_subplot(2,2,[.1 .1],[.1 .05],[.05 .05])

axes(ha(1))
hg1 = line([1369,1369],[-23.97,7.97], 'LineWidth',66,'Color',[0.9,0.9,0.9]);
set(hg1, 'LineWidth',50,'Color',[0.9,0.9,0.9])
hold on
hg2 = line([2470,2470],[-23.97,7.97], 'LineWidth',12,'Color',[0.9,0.9,0.9]);
set(hg2, 'LineWidth',12,'Color',[0.9,0.9,0.9])
hold on
p1 = plot(-mean(CoESfullPOTASi2s,2),'LineWidth',0.8,'Color','#0072BD'),title('(a) CoVaR and CoES series forecasts (Industries \rightarrow System: GAS-CoAS)','position',[1600,-28.5]);
hold on
p2 = plot(-mean(CoVaRfullPOTASi2s,2),'LineWidth',0.8,'Color','r')
hold on
% plot(mean(shenwan,2),'LineWidth',0.3,'Color','#C0C0C0')
p3 = plot(hushen,'LineWidth',0.3,'Color','[0.6,0.6,0.6]')
legend([p1,p2,p3],'CoVaR','CoES','CSI 300','Location','SouthWest','FontSize',8)
set(gca,'XTick',ymlength1);
set(gca,'XTickLabel',dastr1)
set(gca,'YTick',[-24:4:8]);
set(gca,'YTickLabel',{'-24','-20','-16','-12','-8','-4','0','4','8'})
xlim([0 ymlength1(end)])
ylim([-24,8])




axes(ha(2))
hg1 = line([1369,1369],[-23.97,7.97], 'LineWidth',66,'Color',[0.9,0.9,0.9]);
set(hg1, 'LineWidth',50,'Color',[0.9,0.9,0.9])
hold on
hg2 = line([2470,2470],[-23.97,7.97], 'LineWidth',12,'Color',[0.9,0.9,0.9]);
set(hg2, 'LineWidth',12,'Color',[0.9,0.9,0.9])
hold on
p1 = plot(-mean(CoESfullPOTASs2i,2),'LineWidth',0.8,'Color','#0072BD'),title('(b) CoVaR and CoES series forecasts (System \rightarrow Industries: GAS-CoAS)','position',[1600,-28.5]);
hold on
p2 = plot(-mean(CoVaRfullPOTASs2i,2),'LineWidth',0.8,'Color','r')
hold on
p3 = plot(mean(shenwan,2),'LineWidth',0.3,'Color','[0.6,0.6,0.6]')
legend([p1,p2,p3],'CoVaR','CoES','average Industry Index','Location','SouthWest','FontSize',8)
set(gca,'XTick',ymlength1);
set(gca,'XTickLabel',dastr1)
set(gca,'YTick',[-24:4:8]);
set(gca,'YTickLabel',{'-24','-20','-16','-12','-8','-4','0','4','8'})
xlim([0 ymlength1(end)])
ylim([-24,8])



axes(ha(3))
hg1 = line([1369,1369],[-8.97,2.97], 'LineWidth',66,'Color',[0.9,0.9,0.9]);
set(hg1, 'LineWidth',50,'Color',[0.9,0.9,0.9])
hold on
hg2 = line([2470,2470],[-8.97,2.97], 'LineWidth',12,'Color',[0.9,0.9,0.9]);
set(hg2, 'LineWidth',12,'Color',[0.9,0.9,0.9])
hold on
p1 = plot(-mean(CoESfullPOTASi2s-CoVaRfullPOTASi2s,2),'LineWidth',0.8,'Color','[.96,.33,.07'),title('(c) Industries \rightarrow System: GAS-CoAS and CoASc','position',[1600,-10.5]);
hold on
p2 = plot(-mean(CoESfullPOTASi2sc-CoVaRfullPOTASi2sc,2),'LineWidth',0.8,'Color','[.00,.63,.95]')
hold on
p3 = plot(mean(CoESfullPOTASi2s./CoVaRfullPOTASi2s,2),'LineWidth',0.8,'Color','[.49,.18,.56]')
hold on
p4 = plot(mean(CoESfullPOTASi2sc./CoVaRfullPOTASi2sc,2),'LineWidth',0.8,'Color','#77AC30')
hold on
plot(-mean(CoESfullPOTASi2s-CoVaRfullPOTASi2s,2),'LineWidth',0.8,'Color','[.96,.33,.07')
hold on
plot(mean(CoESfullPOTASi2s./CoVaRfullPOTASi2s,2),'LineWidth',0.8,'Color','[.49,.18,.56]')
legend([p1,p2,p3,p4],'CoES-CoVaR(GAS-CoAS)','CoES-CoVaR(CoASc)','CoES/CoVaR(GAS-CoAS)','CoES/CoVaR(CoASc)','Location','SouthWest','FontSize',7.5)
set(gca,'XTick',ymlength1);
set(gca,'XTickLabel',dastr1)
set(gca,'YTick',[-9:3:3]);
set(gca,'YTickLabel',{'-9','-6','-3','0','3'})
xlim([0 ymlength1(end)])
ylim([-9,3])
set(gca,'XGrid','off')



axes(ha(4))
hg1 = line([1369,1369],[-8.97,2.97], 'LineWidth',66,'Color',[0.9,0.9,0.9]);
set(hg1, 'LineWidth',50,'Color',[0.9,0.9,0.9])
hold on
hg2 = line([2470,2470],[-8.97,2.97], 'LineWidth',12,'Color',[0.9,0.9,0.9]);
set(hg2, 'LineWidth',12,'Color',[0.9,0.9,0.9])
hold on
p1 = plot(-mean(CoESfullPOTASs2i-CoVaRfullPOTASs2i,2),'LineWidth',0.8,'Color','[.96,.33,.07'),title('(d) System \rightarrow Industries: GAS-CoAS and CoASc','position',[1600,-10.5]);
hold on
p2 = plot(-mean(CoESfullPOTASs2ic-CoVaRfullPOTASs2ic,2),'LineWidth',0.8,'Color','[.00,.63,.95]')
hold on
p3 = plot(mean(CoESfullPOTASs2i./CoVaRfullPOTASs2i,2),'LineWidth',0.8,'Color','[.49,.18,.56]')
hold on
p4 = plot(mean(CoESfullPOTASs2ic./CoVaRfullPOTASs2ic,2),'LineWidth',0.8,'Color','#77AC30')
hold on
p5 = plot(-mean(CoESfullPOTASs2i-CoVaRfullPOTASs2i,2),'LineWidth',0.8,'Color','[.96,.33,.07');
hold on
p6 = plot(mean(CoESfullPOTASs2i./CoVaRfullPOTASs2i,2),'LineWidth',0.8,'Color','[.49,.18,.56]')
legend([p1,p2,p3,p4],'CoES-CoVaR(GAS-CoAS)','CoES-CoVaR(CoASc)','CoES/CoVaR(GAS-CoAS)','CoES/CoVaR(CoASc)','Location','SouthWest','FontSize',7.5)
set(gca,'XTick',ymlength1);
set(gca,'XTickLabel',dastr1)
set(gca,'YTick',[-9:3:3]);
set(gca,'YTickLabel',{'-9','-6','-3','0','3'})
xlim([0 ymlength1(end)])
ylim([-9,3])

set(gca,'XGrid','off')



