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
figure()
plot(1:60,zeros(60,1),'color','k') % x-axis
hold on
for i = 1:16
    plot(1:60,res50(:,i)/50,'LineStyle',':','color',[1 0 0 0.4])
    hold on
end

% Plot Ind25 errors
for i = 1:13
    plot(1:60,res25(:,i)/25,'LineStyle',':','color',[0 1 0 0.4])
    hold on
end

sh25 = shadedErrorBar(1:60,nanmean(res25,2)/25,nanstd(res25,0,2)/25,{'color','g'},1,1);
sh50 = shadedErrorBar(1:60,mean(res50,2)/50,std(res50,0,2)/50,{'color','r'},1,1);
set(sh25.edge(1),'visible','off')
set(sh25.edge(2),'visible','off')
set(sh50.edge(1),'visible','off')
set(sh50.edge(2),'visible','off')

r25 = plot(1:60,nanmean(res25,2)/25,'color','g','LineWidth',4);
r50 = plot(1:60,mean(res50,2)/50,'color','r','LineWidth',4);
xlabel('time','FontSize',14)
ylabel('error','FontSize',14)
xlim([1 60])
ylim([-1 1])
legend([r50 r25],{'Ind50','Ind25'},'location','best')
title('number evacuated error, individual trials','FontSize',18)
% savefig('figures/evac_error')
% print('figures/evac_error','-dpng','-r300')
%%
figure()
load data/mle_Ind50 rP_hits
h = bar(res50(end,:));
xlabel('trial','fontsize',14)
ylabel('final number evacuated error','fontsize',14)
ybuff=.9;
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
    text(x,y,t,'Color','k','HorizontalAlignment','center')
    ylim([-20 30])
    xlim([0.5 16.5])
end
title('final evacuation error, Ind50 trials', 'fontsize',16)
% savefig('figures/finalerror_50')
% print('figures/finalerror_50','-dpng','-r300')
%%
figure()
load data/mle_Ind25 rP_hits
h = bar(res25(end,:),'m');
xlabel('trial','fontsize',14)
ylabel('final number evacuated error','fontsize',14)
ybuff=1.2;
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
    text(x,y,t,'Color','k','HorizontalAlignment','center')
    ylim([-30 30])
    xlim([0.5 13.5])
end
title('final evacuation error, Ind25 trials', 'fontsize',16)
% savefig('figures/finalerror_25')
% print('figures/finalerror_25','-dpng','-r300')
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


