% This function creates a figure containing subplots of each trial in a
% specified set, where each subplot contains only experimental data, i.e.
% the empirical number evacuated and the observed P_hit trajectory
% Input arguments:
%   z: the numbers of desired trials (leave empty if specific type of trial
%       is desired)
%   saveFig (optional): 1 = saves figures
%   trial_type (optional): string containing the type of trial to be
%       plotted (Ind50, Ind25, Ind5, Grp5, Grp25, FTG5, FTG25, LTG5, LTG25,
%       MV5, MV25, Grp5-S5, Grp5-S25, Grp5-S50, Grp25-S25, Grp25-S50, 
%       FTG5-S5, FTG5-S25, FTG5-S50, FTG25-S25, FTG25-S50, LTG5-S5, LTG5-S25, 
%       LTG5-S50, LTG25-S25, LTG25-S50, MV5-S5, MV5-S25, MV5-S50, MV25-S25, 
%       MV25-S50)
% This function requires the function 'tightfig'

function plot_data(z,saveFig,trial_type)
if nargin < 2
    trial_type = [];
    saveFig = 0;
end

load experiment evapeocumu gameinfo gamestatus missing
gameinfo([30 74 121],39:end) = 1;

if ~isempty(z) % z has priority over trial_type
    trial_type = [];
end

ttype = trial_type;

if ~isempty(trial_type)
    switch trial_type
        case 'Ind50', z = find(gamestatus(:,1)==0&gamestatus(:,2)==50); 
        case 'Ind25', z = find(gamestatus(:,1)==0&gamestatus(:,2)==25);
        case 'Ind5',  z = find(gamestatus(:,1)==0&gamestatus(:,2)==5);
        case 'Grp5',  z = find(gamestatus(:,3)==5); ttype = [];
        case 'Grp25', z = find(gamestatus(:,3)==25); ttype = [];
        case 'FTG5',  z = find(gamestatus(:,1)==1&gamestatus(:,3)==5);
        case 'FTG25', z = find(gamestatus(:,1)==1&gamestatus(:,3)==25);
        case 'LTG5', z = find(gamestatus(:,1)==2&gamestatus(:,3)==5);
        case 'LTG25', z = find(gamestatus(:,1)==2&gamestatus(:,3)==25);
        case 'MV5', z = find(gamestatus(:,1)==3&gamestatus(:,3)==5);
        case 'MV25', z = find(gamestatus(:,1)==3&gamestatus(:,3)==25);
        case 'Grp5-S5', z = find(gamestatus(:,3)==5&gamestatus(:,2)==5); ttype = [];
        case 'Grp5-S25', z = find(gamestatus(:,3)==5&gamestatus(:,2)==25); ttype = [];
        case 'Grp5-S50', z = find(gamestatus(:,3)==5&gamestatus(:,2)==50); ttype = [];
        case 'Grp25-S25', z = find(gamestatus(:,3)==25&gamestatus(:,2)==25); ttype = [];
        case 'Grp25-S50', z = find(gamestatus(:,3)==25&gamestatus(:,2)==50); ttype = [];
        case 'FTG5-S5', z = find(gamestatus(:,1)==1&gamestatus(:,2)==5&gamestatus(:,3)==5); ttype = 'FTG5';
        case 'FTG5-S25', z = find(gamestatus(:,1)==1&gamestatus(:,2)==25&gamestatus(:,3)==5); ttype = 'FTG5';
        case 'FTG5-S50', z = find(gamestatus(:,1)==1&gamestatus(:,2)==50&gamestatus(:,3)==5); ttype = 'FTG5';
        case 'FTG25-S25', z = find(gamestatus(:,1)==1&gamestatus(:,2)==25&gamestatus(:,3)==25); ttype = 'FTG25';
        case 'FTG25-S50', z = find(gamestatus(:,1)==1&gamestatus(:,2)==50&gamestatus(:,3)==25); ttype = 'FTG25';
        case 'LTG5-S5', z = find(gamestatus(:,1)==2&gamestatus(:,2)==5&gamestatus(:,3)==5); ttype = 'LTG5';
        case 'LTG5-S25', z = find(gamestatus(:,1)==2&gamestatus(:,2)==25&gamestatus(:,3)==5); ttype = 'LTG5';
        case 'LTG5-S50', z = find(gamestatus(:,1)==2&gamestatus(:,2)==50&gamestatus(:,3)==5); ttype = 'LTG5';
        case 'LTG25-S25', z = find(gamestatus(:,1)==2&gamestatus(:,2)==25&gamestatus(:,3)==25); ttype = 'LTG25';
        case 'LTG25-S50', z = find(gamestatus(:,1)==2&gamestatus(:,2)==50&gamestatus(:,3)==25); ttype = 'LTG25';
        case 'MV5-S5', z = find(gamestatus(:,1)==3&gamestatus(:,2)==5&gamestatus(:,3)==5); ttype = 'MV5';
        case 'MV5-S25', z = find(gamestatus(:,1)==3&gamestatus(:,2)==25&gamestatus(:,3)==5); ttype = 'MV5';
        case 'MV5-S50', z = find(gamestatus(:,1)==3&gamestatus(:,2)==50&gamestatus(:,3)==5); ttype = 'MV5';
        case 'MV25-S25', z = find(gamestatus(:,1)==3&gamestatus(:,2)==25&gamestatus(:,3)==25); ttype = 'MV25';
        case 'MV25-S50', z = find(gamestatus(:,1)==3&gamestatus(:,2)==50&gamestatus(:,3)==25); ttype = 'MV25';
    end
end

z = z(~ismember(z,missing));
rP_hits = gameinfo(z,:);
rP_hits = 10.^floor(log10(abs(rP_hits)));
rP_hits = floor(gameinfo(z,:)./rP_hits).*rP_hits;
rP_hits(isnan(rP_hits)) = 0;
rP_hits = round(rP_hits,1);

figure()
for i = 1:length(z)
    if isempty(trial_type)
        if gamestatus(z(i),1)==0
            ttype = 'Ind';
            ttype = strcat(ttype,num2str(gamestatus(z(i),2)));
        elseif gamestatus(z(i),1)==1
            ttype = 'FTG';
            ttype = strcat(ttype,num2str(gamestatus(z(i),3)));
        elseif gamestatus(z(i),1)==2
            ttype = 'LTG';
            ttype = strcat(ttype,num2str(gamestatus(z(i),3)));
        elseif gamestatus(z(i),1)==3
            ttype = 'MV';
            ttype = strcat(ttype,num2str(gamestatus(z(i),3)));
        end         
    end
        
    subplot(ceil(length(z)/4),4,i)
    hold on
    area(1:60,evapeocumu(z(i),:),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')    
    [ax,shelter,phits]=plotyy(1:60,ones(60,1)*gamestatus(z(i),2),1:60,rP_hits(i,:));
    set(phits,'color',[0 1 0]); set(phits,'LineWidth',1); set(phits,'LineStyle','-');
    set(shelter,'color',[0 0 1]); set(shelter,'LineWidth',4); set(shelter,'LineStyle',':');
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).XLabel.String = 'time';
    ax(1).YTick = 0:10:50; ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50]; ax(2).YLim = [0 1];
    ax(1).XLim = [1 60]; ax(2).XLim = [1 60];
    ax(1).FontSize = 11; ax(2).FontSize = 11;
    title(['Trial ' num2str(z(i)) ' (' ttype ')'],'FontSize',14)
end


set(gcf,'units','pixels')
set(gcf,'position',[0 0 1300 800])
tightfig;
pos=get(gcf,'position');
set(gcf,'position',[pos(1:3) 800])
tightfig;

if saveFig
    savefig(['figures/' trial_type])
    print(['figures/' trial_type],'-dpdf','-r300')
    print(['figures/' trial_type],'-dsvg','-r300')
end


end