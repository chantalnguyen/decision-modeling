% Plot H (P_hits observed), J (evacuations), and H/J (evacuation rate) for
% Ind50 trials

% desired trial #s are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144
z = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
space = 50;

% misses 28 29 64 67 113 125 130 143 144
% hits 19 49 76 102 108 112 123

% import data 
load evacuate_data

% check to see if any bad trials are in z
try assert(~any(ismember(z,missing)))
catch err
    error('Missing trials requested!');
end

% round rP_hits down to nearest tenth, since participants only see discrete
% values of P_hit
round_P_hits = 10.^floor(log10(abs(gameinfo)));
round_P_hits = floor(gameinfo./round_P_hits).*round_P_hits;
round_P_hits(isnan(round_P_hits))=0;
round_P_hits = round(round_P_hits,1);

evacuateTime(evacuateTime == 0) = 1; % evacuating at 0 = evacuating at 1

rP_hits = zeros(60,length(z)); % rounded P_hit trajectories for each trial
evac = rP_hits; % empirical cumulative evacuations for each trial
evacTimes = zeros(50,length(z)); % times of evacuation decisions for each trial
evacPhits = evacTimes; % P_hits of evacuation decisions for each trial

% this is our data set
for i = 1:length(z)
    rP_hits(:,i) = round_P_hits(z(i),:);
    evac(:,i) = evapeocumu(z(i),:); 
    evacTimes(:,i) = evacuateTime(:,z(i));
    evacPhits(:,i) = evacuateProb(:,z(i));
end

% round evacuateProb values down to nearest tenth
evacPhits2 = 10.^floor(log10(abs(evacPhits)));
evacPhits2 = floor(evacPhits./evacPhits2).*evacPhits2;
evacPhits2(isnan(evacPhits2))=0;
evacPhits = round(evacPhits2,1);
clear evacPhits2;

P_hit_range=0:0.1:1;

% Number of times participants at home observe P_hit values
H = zeros(length(z),length(P_hit_range));

% Number of times participants decide to evacuate at P_hit values
J = zeros(length(z),length(P_hit_range));

bins = -0.05:0.1:1.05;


% for j = 1:length(z) % iterate through each trial
%     for i = 1:50 % iterate through each individual
%         indvEvacTime = evacTimes(i,j); % this is the time the individual (i) evacuated in trial (j)
%         if indvEvacTime == -1 % no decision; ALL rP_hits seen until end of trial (P_hit = 1 or 0)
%             assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
%             h1 = histcounts(rP_hits(1:find(gameinfo(z(j),:)==gameinfo(z(j),end),1,'first'),j),bins);
%         else % decision; count P_hits seen until time of decision
%             assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
%             h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
%         end
%         H(j,:) = H(j,:)+h1;
%     end
% end

if space == 50
    for j = 1:length(z) % iterate through each trial
        for i = 1:50 % iterate through each individual
            indvEvacTime = evacTimes(i,j); % this is the time the individual (i) evacuated in trial (j)
            if indvEvacTime == -1 % no decision; ALL rP_hits seen until end of trial (P_hit = 1 or 0)
                assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
                h1 = histcounts(rP_hits(1:find(gameinfo(z(j),:)==gameinfo(z(j),end),1,'first'),j),bins);
            else % decision; count P_hits seen until time of decision
                assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
                h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
            end
            H(j,:) = H(j,:)+h1;
        end
    end
else
    for j = 1:length(z)
        for i = 1:50
            indvEvacTime = evacTimes(i,j);
            if indvEvacTime == -1
                if evac(end,j) < space
                    h1 = histcounts(rP_hits(1:find(gameinfo(z(j),:)==gameinfo(z(j),end),1,'first'),j),bins);
                else
                    h1 = histcounts(rP_hits(1:find(evac(:,j)==space,1,'first'),j),bins);
                end
            else
                h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
            end
            H(j,:) = H(j,:) + h1;
        end
    end
end


for i = 1:length(z) % iterate over all trials
    J(i,:) = J(i,:) + histcounts(evacPhits(:,i),bins);
end

Theta = J'./H';
%Theta = (J'+1)./(H'+2);
%% Plot H
figure()
h = pcolor(-0.1:0.1:1,1:length(z)+1,padarray(H,[1 1],mean(mean(H)),'post'));
set(h,'EdgeColor','none');
colormap(parula(max(max(H))-min(min(H))))
set(gca,'xtick',-0.1+0.05:0.1:1,'xticklabel',(-0.1:0.1:1)+0.1);
set(gca,'ytick',2.5:2:length(z)+1,'yticklabel',2:2:length(z));
set(gca,'fontsize',16);
b=colorbar;
b.Label.String = 'Times observed';
b.Label.FontSize = 18;
xlabel('Disaster likelihood (P_{hit})','fontsize',18)
ylabel('Trial','fontsize',18)
title('Times P_{hit} value observed H, Ind50 trials','fontsize',18)
savefig('figures/Ind50_H')
print('figures/Ind50_H','-dsvg','-r300')
%% Plot J
figure()
h = pcolor(-0.1:0.1:1,1:length(z)+1,padarray(J,[1 1],mean(mean(J)),'post'));
set(h,'EdgeColor','none');
colormap(parula(max(max(J))+1-min(min(J))))
set(gca,'xtick',-0.1+0.05:0.1:1,'xticklabel',(-0.1:0.1:1)+0.1);
set(gca,'ytick',2.5:2:length(z)+1,'yticklabel',2:2:length(z));
set(gca,'fontsize',16);
b=colorbar;
b.Label.String = 'Evacuations';
b.Label.FontSize = 18;
b.Ticks = [0 5.25 10.125 15 19.75 24.625 30];
b.TickLabels={'0','5','10','15','20','25','30'};
xlabel('Disaster likelihood (P_{hit})','fontsize',18)
ylabel('Trial','fontsize',18)
title('Empirical evacuations J, Ind50 trials','fontsize',18)
savefig('figures/Ind50_J')
print('figures/Ind50_J','-dsvg','-r300')
%% Plot J/H
figure()
h=bar(0:0.1:1,Theta,'stacked');
set(h,'edgecolor','k');
xlim([-0.05 1.05])
set(gca,'fontsize',16)
set(gca,'xtick',0:0.1:1)
set(gca,'ticklength',[0 0])
xlabel('Disaster likelihood (P_{hit})','fontsize',18)
ylabel('Evacuation rate','fontsize',18)
title('Empirical evacuation rate J/H, Ind50 trials', 'fontsize',18)
colormap(jet(length(z)));
legs = cell(length(z),1);
for i = 1:length(z)
    legs{i}=['Ind50-' num2str(i)];
end
legend(legs,'location','northwest','fontsize',14)

%% Plot J/H, all one color
figure()
h=bar(0:0.1:1,nansum(Theta'),'b');
% set(h,'edgecolor','none');
xlim([-0.05 1.05])
set(gca,'fontsize',16)
set(gca,'xtick',0:0.1:1)
set(gca,'ticklength',[0 0])
xlabel('Disaster likelihood (P_{hit})','fontsize',18)
ylabel('Evacuation rate','fontsize',18)
title('Empirical evacuation rate', 'fontsize',18)
% colormap(jet(length(z)));
% legs = cell(length(z),1);
% for i = 1:length(z)
%     legs{i}=['Ind50-' num2str(i)];
% end
% legend(legs,'location','northwest','fontsize',14)
% savefig('figures/Ind50_JH')
% print('figures/Ind50_JH','-dsvg','-r300')
%% Plot cumulative J/H
cTheta = nanmean(Theta,2);
cTheta = cumsum(cTheta)/max(cumsum(cTheta));
cJ = sum(J);
cJ = cumsum(cJ)/max(cumsum(cJ));
Theta2 = (J+1)'./(H+2)'; % actual mean
cTheta2 = nanmean(Theta2,2);
cTheta2 = cumsum(cTheta2)/max(cumsum(cTheta2));
figure()
plot(0:0.1:1,cTheta,'linewidth',2); hold on;
% plot(0:0.1:1,cJ,'linewidth',2) 
plot(0:0.1:1,0.8824*(0:0.1:1).^5.4053,'linewidth',2)
plot(0:0.1:1,0.0394*(0:0.1:1).^0.8062,'linewidth',2)
% plot(0:0.1:1,cTheta2,'linestyle','--')
set(gca,'fontsize',16)
xlabel('Disaster likelihood (P_{hit})','fontsize',18)
% ylabel('Cumulative evacuation rate','fontsize',18)
title('Cumulative evacuation rate, Ind50 trials','fontsize',18)
legend({'Cumulative evacuation rate','Decision model'},'Location','NorthWest')

