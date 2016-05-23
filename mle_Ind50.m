% Determine parameters of power-law decision model for Ind50 games with
% maximum likelihood estimation, and solve master equation for each trial

% desired trial #s are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144
z = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];

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
H = zeros(50,length(P_hit_range));

% Number of times participants decide to evacuate at P_hit values
J = zeros(50,length(P_hit_range));

bins = -0.05:0.1:1.05;

for i = 1:50 % iterate through each individual
    for j = 1:length(z) % iterate through each trial
        indvEvacTime = evacTimes(i,j); % this is the time the individual (j) evacuated in trial (r)
        if indvEvacTime == -1 % no decision; ALL rP_hits seen until end of trial (P_hit = 1 or 0)
            assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
            h1 = histcounts(rP_hits(1:find(rP_hits(:,j)==rP_hits(end,j),1,'first'),j),bins);
        else % decision; count P_hits seen until time of decision
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
powerfun = @(theta_p)llfun(H(2:end),J(2:end),theta_p,P_hit_range(2:end));

options = optimoptions(@fminunc,'MaxFunEvals',10000);

% minimize (negative of) log-likelihood function
[theta,~] = fminunc(powerfun,theta_0,options);

power_model = @(Ph,theta_p) theta_p(1)*Ph.^theta_p(2);

% determine standard deviations from boostrapping
[bootstat,bootsam]=bootstrp(1000,@(z)mle_fit(z,gameinfo,evapeocumu,evacuateTime,evacuateProb),z);
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
    endTimes(i) = find(rP_hits(:,i)==rP_hits(end,i),1,'first');
end

P_Ind50 = zeros(60,16);
Probs_Ind50 = cell(16,1);
mProbs_Ind50 = cell(16,1);
stdevs_Ind50 = cell(16,1);
T_Ind50 = cell(16,1);
for i = 1:length(z)
    [Ptest,Ttest,PPtest] = mastereq(power_model(rP_hits(:,i),theta),endTimes(i));
    mProbs_Ind50{i} = Ptest;
    T_Ind50{i} = Ttest;
    temp = interp1(Ttest,Ptest,1:1:endTimes(i));
    temp(endTimes(i)+1:60)=temp(end);
    P_Ind50(:,i) = temp;
    Probs_Ind50{i} = PPtest;
    stdevs_Ind50{i} = zeros(size(Ttest));
    for j = 1:length(Ttest)
        stdevs_Ind50{i}(j) = sqrt(((0:50)-Ptest(j)).^2*(PPtest(j,:)'));
    end
end

% calculate root mean square error
rss_Ind50 = zeros(length(z),1);
for i=1:length(z)
    for j = 1:endTimes(i)
        rss_Ind50(i)=rss_Ind50(i)+(P_Ind50(j,i)-evac(j,i))^2;
    end
    rss_Ind50(i)=rss_Ind50(i)/endTimes(i);
end
mse_Ind50 = mean(rss_Ind50);
rmse_Ind50 = sqrt(mse_Ind50);

% save workspace
clear bins temp i j
save('data/mle_Ind50.mat')
    