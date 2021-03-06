% This function creates a figure containing subplots of each trial in a
% specified set, where each subplot contains the empirical number
% evacuated, the average of the modeled probability distribution, the
% 3-sigma confidence interval for the model average, the cross-validation
% prediction, and the empirical P_hit trajectory
% The lower left subplot is enlarged to include axis labels and a legend
% Input arguments:
%   input_data: a string containing the name of the basic input dataset
%       with model prediction
%   input_cv: a string containing the name of the corresponding dataset
%       with cross-validation prediction
%   trial_type: a string containing the type of trial being plotted, e.g.
%       Ind50, Ind25, LTG50
%   saveFig (optional): 1 = saves figures
% This function requires the functions 'tightfig' and 'shadedErrorBar'
% This really only works if there are between 10 and 17 trials

function plot_opt_trials(input_data,input_bayes,input_stat,trial_type,saveFig)
if nargin < 5
    saveFig = 0;
end


load(input_data);

temp = who('T_*');      Times = eval(temp{1});
temp = who('mProbs_Ind50'); meanProbs = eval(temp{1});

load(input_bayes);
load(input_stat);
temp = who('P_bayes*'); P_bayes = eval(temp{1});
temp = who('P_stat*');  P_stat = eval(temp{1});
temp = who('rP_hits');  rP_hits = eval(temp{1});
temp = who('z');        z = eval(temp{1});
temp = who('endTimes'); endTimes = eval(temp{1});

if strncmpi(trial_type,'LTG',3)
    load evacuate_data evacuateTime
    evac = zeros(60,length(z));
    for i = 1:length(z)
        decTimes = evacuateTime(:,z(i));
        decTimes(decTimes==-1) = NaN;
        decTimes = sort(decTimes);
        N = histcounts(decTimes,0.5:1:60.5);
        evac(:,i) = cumsum(N);
    end
else
    temp = who('evac');     evac = eval(temp{1});
end

if length(z) >=10 && length(z) <=13
    numrows = 4; firstPlots = 1:8; lastPlot = [11 12 15 16];
else
    numrows = 5; firstPlots = 1:12; lastPlot = [15 16 19 20];
end

figure()
% plot first 2 or 3 rows
for i = firstPlots
    subplot(numrows,4,i) 
    hold on
    area(1:60,evac(:,i),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')    
    plot(Times{i},meanProbs{i},'color','k','LineWidth',5,'LineStyle',':');
    plot(1:60,P_stat(:,i),'Linewidth',5,'color',[1 .5 0]);
    [ax,bayes,phits]=plotyy(1:60,P_bayes(:,i),1:60,rP_hits(:,i));
    set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
    set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(i)]; ax(2).XLim = [1 endTimes(i)];
    ax(1).FontSize = 11; ax(2).FontSize = 11;
    title(['Trial ' trial_type '-' num2str(i)],'FontSize',14)
end

% plot bottom right plot
lP = length(z);
subplot(numrows,4,lastPlot)
    hold on
    evac_obs = area(1:60,evac(:,lP),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')
    model = plot(Times{i},meanProbs{i},'color','k','LineWidth',5,'LineStyle',':'); 
    statopt = plot(1:60,P_stat(:,i),'Linewidth',5,'color',[1 .5 0]);
    [ax,bayes,phits]=plotyy(1:60,P_bayes(:,lP),1:60,rP_hits(:,lP));
    set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
    set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(lP)]; ax(2).XLim = [1 endTimes(lP)];
    ax(1).FontSize = 13; ax(2).FontSize = 13;
    ax(2).YLabel.FontSize=18; ax(2).YLabel.String = 'Likelihood';
    xlabel('Time step','FontSize',18)
    ylabel('Evacuations','FontSize',18)
    title(['Trial ' trial_type '-' num2str(lP)],'FontSize',14)
    legend([evac_obs, model,bayes,statopt,phits],{'   Data','   Model average','   Bayesian optimal','   Static optimal','   P_{hit}'},'location','northwest','fontsize',17)
   
% plot remaining plots
if length(z) == 10 || length(z) == 14
    plotInd = firstPlots(end)+2;
    subplot(numrows,4,plotInd);
    hold on
    area(1:60,evac(:,plotInd),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')    
    plot(Times{plotInd},meanProbs{plotInd},'color','k','LineWidth',5,'LineStyle',':');
    plot(1:60,P_stat(:,plotInd),'Linewidth',5,'color',[1 .5 0]);    
    [ax,bayes,phits]=plotyy(1:60,P_bayes(:,plotInd),1:60,rP_hits(:,plotInd));
    set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
    set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(plotInd)]; ax(2).XLim = [1 endTimes(plotInd)];
    ax(1).FontSize = 11; ax(2).FontSize = 11;
    title(['Trial ' trial_type '-' num2str(plotInd)],'FontSize',14)
else
    for j = (firstPlots(end)+1):(firstPlots(end)+2)
        subplot(numrows,4,j);
        hold on
        area(1:60,evac(:,j),'FaceColor',[192/255 192/255 192/255]);
        set(gca,'Layer','top')    
        plot(Times{j},meanProbs{j},'color','k','LineWidth',5,'LineStyle',':');
        plot(1:60,P_stat(:,j),'Linewidth',5,'color',[1 .5 0]);
        [ax,bayes,phits]=plotyy(1:60,P_bayes(:,j),1:60,rP_hits(:,j));
        set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
        set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
        set(ax(1),'box','on')
        set(ax(2),'YColor','g')
        ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
        ax(1).YColor = 'k';
        ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
        ax(1).XLim = [1 endTimes(j)]; ax(2).XLim = [1 endTimes(j)];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        title(['Trial ' trial_type '-' num2str(j)],'FontSize',14)
    end
    if length(z) == 12 || length(z) == 16
        plotInd = length(z) - 1;
        subplot(numrows,4,plotInd+3);
        hold on
        area(1:60,evac(:,plotInd),'FaceColor',[192/255 192/255 192/255]);
        set(gca,'Layer','top')        
        plot(Times{plotInd},meanProbs{plotInd},'color','k','LineWidth',5,'LineStyle',':');
        plot(1:60,P_stat(:,plotInd),'Linewidth',5,'color',[1 .5 0]);        
        [ax,bayes,phits]=plotyy(1:60,P_bayes(:,plotInd),1:60,rP_hits(:,plotInd));
        set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
        set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
        set(ax(1),'box','on')
        set(ax(2),'YColor','g')
        ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
        ax(1).YColor = 'k';
        ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
        ax(1).XLim = [1 endTimes(plotInd)]; ax(2).XLim = [1 endTimes(plotInd)];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        title(['Trial ' trial_type '-' num2str(plotInd)],'FontSize',14)
    elseif length(z) == 13 || length(z) == 17
        for k = (length(z)-2):(length(z)-1)
            subplot(numrows,4,k+2);
            hold on
            area(1:60,evac(:,k),'FaceColor',[192/255 192/255 192/255]);
            set(gca,'Layer','top')        
            plot(Times{i},meanProbs{i},'color','k','LineWidth',5,'LineStyle',':');
            plot(1:60,P_stat(:,i),'Linewidth',5,'color',[1 .5 0]);
            [ax,bayes,phits]=plotyy(1:60,P_bayes(:,k),1:60,rP_hits(:,k));
            set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
            set(bayes,'color',[0 1 1 0.4]); set(bayes,'LineWidth',5);
            set(ax(1),'box','on')
            set(ax(2),'YColor','g')
            ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
            ax(1).YColor = 'k';
            ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
            ax(1).XLim = [1 endTimes(k)]; ax(2).XLim = [1 endTimes(k)];
            ax(1).FontSize = 11; ax(2).FontSize = 11;
            title(['Trial ' trial_type '-' num2str(k)],'FontSize',14)
        end
    end
end


set(gcf,'units','pixels')
set(gcf,'position',[0 0 1300 800])
tightfig;
pos=get(gcf,'position');
set(gcf,'position',[pos(1:3) 800])
tightfig;
% set(gcf,'units','centimeters')
% set(gcf,'position',[pos(1:2) 40 28])
% set(gcf,'papersize',[40 28])

if saveFig
    savefig(['figures/opt_' trial_type])
    print(['figures/opt_' trial_type],'-dpdf','-r300')
    print(['figures/opt_' trial_type],'-dsvg','-r300')
end


end