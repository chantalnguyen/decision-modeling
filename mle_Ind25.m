  % Determine parameters of power-law decision model for Ind25 games with
% maximum likelihood estimation, and solve master equation for each trial
% To calculate H, count P_hits seen until shelter is full
% Decision model is dependent on P_hit and initial shelter capacity (25)

% desired trial #s
z = [34 36 52 61 65 66 71 82 95 139 148 152 160];
space = 25;

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
H = zeros(50,length(P_hit_range));

% Number of times participants decide to evacuate at P_hit values
J = zeros(50,length(P_hit_range));

bins = -0.05:0.1:1.05;

for i = 1:50 % iterate through each individual
    for j = 1:length(z) % iterate through each trial
        indvEvacTime = evacTimes(i,j); % this is the time the individual (i) evacuated in trial (j)
        if indvEvacTime == -1 % no decision; ALL rP_hits seen until when shelter is full NOT end of trial (P_hit = 1 or 0) 
                assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
            if evac(end,j) < space % if shelter did not fill, count all rP_hits seen
                h1 = histcounts(rP_hits(1:find(gameinfo(z(j),:)==gameinfo(z(j),end),1,'first'),j),bins);
            else % if shelter fills up, count rP_hits seen until shelter is full
                 % honestly, this would cause a discrepancy if
                 % participants decide when the shelter is already
                 % filled
                h1 = histcounts(rP_hits(1:find(evac(:,j)==space,1,'first'),j),bins); 
            end
        else 
                assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
            h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
        end
        H(i,:) = H(i,:)+h1;
    end
end


for i = 1:50 % iterate over all individials
    J(i,:) = J(i,:) + histcounts(evacPhits(i,:),bins);
end


H = sum(H)';
J = sum(J)';

theta_0 = [0.8; 10];
powerfun = @(theta_i)llfun(H(2:end),J(2:end),theta_i,P_hit_range(2:end),space);

options = optimoptions(@fminunc,'MaxFunEvals',10000);

% minimize (negative of) log-likelihood function
[theta,~] = fminunc(powerfun,theta_0,options);

power25model = @(Ph,space,theta_i) theta_i(1)*Ph.^(theta_i(2)/space);

% determine standard deviations from boostrapping
[bootstat,bootsam]=bootstrp(1000,@(z)mle_fit(z,gameinfo,evapeocumu,evacuateTime,evacuateProb,space),z);
std_theta = std(bootstat);

% plot decision model
% figure()
% Phits = 0:0.01:1;
% plot(Phits,power_model(Phits,theta_power),'LineWidth',3);
% xlabel('P_{hit}','FontSize',16)
% ylabel('probability of evacuation q(P_{hit})','FontSize',16)
% title('decision model for individual games','FontSize',18)
% set(gca,'FontSize',12);

%% solve master equation

% calculate endtimes of each trial
endTimes = zeros(size(rP_hits,2),1);
for i = 1:size(rP_hits,2)
    endTimes(i) = find(gameinfo(z(i),:)==gameinfo(z(i),end),1,'first'); % use unrounded values
end


P_Ind25 = zeros(60,length(z)); % mean cumulative evacuations for each trial
Probs_Ind25 = cell(length(z),1); % full probability distributions for each trial
mProbs_Ind25 = cell(length(z),1); % (raw) mean cumulative evacuations for each trial
stdevs_Ind25 = cell(length(z),1); % standard deviation of mean cumulative evacuations
T_Ind25 = cell(length(z),1); % ODE time steps for each trial
for i = 1:length(z)
    [Ptest,Ttest,PPtest] = mastereq(power25model(rP_hits(:,i),space,theta),endTimes(i),space);
    mProbs_Ind25{i} = Ptest;
    T_Ind25{i} = Ttest;
    temp = interp1(Ttest,Ptest,1:1:endTimes(i)); 
    temp(endTimes(i)+1:60)=temp(end);
    P_Ind25(:,i) = temp;
    Probs_Ind25{i} = PPtest;
    stdevs_Ind25{i} = zeros(size(Ttest));
    for j = 1:length(Ttest)
        stdevs_Ind25{i}(j) = sqrt(((0:50)-Ptest(j)).^2*(PPtest(j,:)'));
    end
end

% calculate root mean square error
rss_Ind25 = zeros(length(z),1);
for i=1:length(z)
    for j = 1:endTimes(i)
        rss_Ind25(i)=rss_Ind25(i)+(P_Ind25(j,i)-evac(j,i))^2;
    end
    rss_Ind25(i)=rss_Ind25(i)/endTimes(i);
end
mse_Ind25 = mean(rss_Ind25);
rmse_Ind25 = sqrt(mse_Ind25);

% save workspace
clear bins temp i j
save('data/mle_Ind25.mat')

%% fit trials sequentially in groups
% thetas = zeros(7,2);
% z = [34 36 52 61 65 66 71 82 95 139 148 152 160];
% space = 25;
% load evacuate_data
% for i = 1:length(thetas)
%     theta = mle_fit(z(i:i+5),gameinfo,evapeocumu,evacuateTime,evacuateProb,space);
%     thetas(i,:) = theta;
% end
% %%
% power_model = @(Ph,space,theta_i) theta_i(1)*Ph.^(theta_i(2)/space);
% % plot decision model
% colors = flipud(jet(length(thetas)));
% legs = cell(length(thetas),1);
% figure('position',[0 0 375 281.25])
% for i = 1:length(thetas)
%     Phits = 0:0.01:1;
%     plot(Phits,power_model(Phits,space,thetas(i,:)),'LineWidth',2,'color',colors(i,:));
%     hold on
%     xlabel('Disaster likelihood, P_{hit}','FontSize',12)
%     ylabel('Probability of evacuation, q(P_{hit})','FontSize',12)
%     title('Decision model for Ind25 trials','FontSize',12)
%     set(gca,'FontSize',12);
%     set(gca,'xtick',0:0.1:1);
%     legs{i} = ['  a = ' num2str(thetas(i,1)) ', b = ' num2str(thetas(i,2))];
% end
% leg=legend(legs,'location','northwest');
% set(leg,'fontsize',12);

    