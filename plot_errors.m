% Calculate and plot absolute evacuation errors (# evacuated in model - #
% evacuated in experiment) for Ind50 and Ind25 trials

load data/mle_Ind50 P_Ind50 evac
evac50 = evac;
load data/mle_Ind25 P_Ind25 evac
evac25 = evac;
clear evac;

res50 = P_Ind50 - evac50;
res25 = P_Ind25 - evac25;

% Plot Ind50 errors
figure('position',[0 0 1200 800])
subplot(2,2,1)
plot(1:60,zeros(60,1),'color','k') % x-axis
hold on
for i = 1:16
    plot(1:60,res50(:,i)/50,'LineStyle',':','color',[1 0 0 0.4])
    hold on
end
sh50 = shadedErrorBar(1:60,mean(res50,2)/50,std(res50,0,2)/50,{'color','r'},1,1);
set(sh50.edge(1),'visible','off'); set(sh50.edge(2),'visible','off');
r50 = plot(1:60,mean(res50,2)/50,'color','r','LineWidth',4);
xlabel('time','FontSize',14)
ylabel('error','FontSize',14)
set(gca,'fontsize',14)
xlim([1 60])
ylim([-1 1])
title('number evacuated error, Ind50 trials','FontSize',14)

% Plot Ind25 errors
subplot(2,2,2)
plot(1:60,zeros(60,1),'color','k') % x-axis
hold on
for i = 1:13
    plot(1:60,res25(:,i)/25,'LineStyle',':','color',[0 0.8 0 0.6])
    hold on
end
sh25 = shadedErrorBar(1:60,nanmean(res25,2)/25,nanstd(res25,0,2)/25,{'color','g'},1,1);
set(sh25.edge(1),'visible','off'); set(sh25.edge(2),'visible','off');
r25 = plot(1:60,nanmean(res25,2)/25,'color','g','LineWidth',4);
xlabel('time','FontSize',14)
ylabel('error','FontSize',14)
set(gca,'fontsize',14)
xlim([1 60])
ylim([-1 1])
title('number evacuated error, Ind25 trials','FontSize',14)

subplot(2,2,3)
load data/mle_Ind50 rP_hits
bar(res50(end,:));
xlabel('trial','fontsize',14)
ylabel('final number evacuated error','fontsize',14)
ybuff=1.5;
XDATA=get(get(gca,'Children'),'XData');
YDATA=get(get(gca,'Children'),'YData');
for j=1:size(XDATA,2)
    x=XDATA(j);
    if YDATA(j)>0
        y=YDATA(j)+ybuff;
    else
        y = YDATA(j)-ybuff;
    end
    if rP_hits(end,j) == 1
        t = 'hit';
    else
        t = 'miss';
    end
    text(x,y,t,'Color','k','HorizontalAlignment','center','fontsize',11)
end
xlim([0.25 16.75])
ylim([-20 30])
set(gca,'xtick',2:2:16)
set(gca,'fontsize',14)
title('final evacuation error, Ind50 trials', 'fontsize',14)

subplot(2,2,4)
load data/mle_Ind25 rP_hits
bar(res25(end,:),'m');
xlabel('trial','fontsize',14)
ylabel('final number evacuated error','fontsize',14)
ybuff=1.5;
XDATA=get(get(gca,'Children'),'XData');
YDATA=get(get(gca,'Children'),'YData');
for j=1:size(XDATA,2)
    x=XDATA(j);
    if YDATA(j)>0
        y=YDATA(j)+ybuff;
    else
        y = YDATA(j)-ybuff;
    end
    if rP_hits(end,j) == 1
        t = 'hit';
    else
        t = 'miss';
    end
    text(x,y,t,'Color','k','HorizontalAlignment','center','fontsize',11)
end
xlim([0.25 13.75])
ylim([-30 30])
set(gca,'fontsize',14)
title('final evacuation error, Ind25 trials', 'fontsize',14)
tightfig;
%%
load data/mle_Ind50 rP_hits
figure()
final_error_hit = mean(res50(end,rP_hits(end,:)==1));
final_error_miss = mean(res50(end,rP_hits(end,:)==0));
h2 = bar([1 2],[final_error_hit final_error_miss],'r');
hold on
errorbar([final_error_hit final_error_miss],[std(res50(end,rP_hits(end,:)==1)),std(res50(end,rP_hits(end,:)==0))],'linestyle','none','color','k')
xlim([0.5 2.5])
ylim([-10 15])
set(gca,'XTickLabel',{'hits','misses'},'fontsize',40)
ylabel('mean final evacuation error','fontsize',37)
tightfig;
set(gcf,'position',[0 0 535 490])
tightfig;
% savefig('figures/meanfinalerror_50')
% print('figures/meanfinalerror_50','-dpng','-r300')
%%
load data/mle_Ind25 rP_hits
figure()
final_error_hit = mean(res25(end,rP_hits(end,:)==1));
final_error_miss = mean(res25(end,rP_hits(end,:)==0));
bar([1 2],[final_error_hit final_error_miss],'r');
hold on
errorbar([final_error_hit final_error_miss],[std(res25(end,rP_hits(end,:)==1)),std(res25(end,rP_hits(end,:)==0))],'linestyle','none','color','k')
xlim([0.5 2.5])
ylim([-10 15])
set(gca,'XTickLabel',{'hits','misses'},'fontsize',40) 
ylabel('mean final evacuation error','fontsize',37) 
tightfig;
set(gcf,'position',[0 0 535 490])
tightfig;
% savefig('figures/meanfinalerror_25')
% print('figures/meanfinalerror_25','-dpng','-r300')


