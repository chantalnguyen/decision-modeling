% Perform cross-validation for Ind50 games and solve master equation for each trial

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
evacTimes = zeros(50,length(z)); % times of evacuation DECISIONS for each trial
evacPhits = evacTimes; % P_hits of evacuation DECISIONS for each trail

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

eTimes = zeros(50,length(z)-1,length(z));
ePhits = eTimes;
looPhits = zeros(60,length(z)-1,length(z));

eTimes(:,:,1) = evacTimes(:,2:end);
ePhits(:,:,1) = evacPhits(:,2:end);
eTimes(:,:,end) = evacTimes(:,1:end-1);
ePhits(:,:,end) = evacPhits(:,1:end-1);
looPhits(:,:,1) = rP_hits(:,2:end);
looPhits(:,:,end) = rP_hits(:,1:end-1);

for i = 2:length(z)
    eTimes(:,:,i) = horzcat(evacTimes(:,1:i-1),evacTimes(:,i+1:end));
    ePhits(:,:,i) = horzcat(evacPhits(:,1:i-1),evacPhits(:,i+1:end));
    looPhits(:,:,i) = horzcat(rP_hits(:,1:i-1),rP_hits(:,i+1:end));
end

P_hit_range=0:0.1:1;
P_hit_range=P_hit_range';

% Number of times participants at home observe P_hit values
H = zeros(50,length(P_hit_range),length(z));
Htrue = zeros(50,length(P_hit_range));

% Number of times participants evacuate at P_hit values
J = zeros(50,length(P_hit_range),length(z));
Jtrue = zeros(50,length(P_hit_range));

bins = -0.05:0.1:1.05;

for n = 1:length(z)
    for i = 1:50 % iterate through each individual
        for j = 1:length(z)-1 % iterate through each trial
            indvEvacTime = eTimes(i,j,n); % this is the time the individual (j) evacuated in trial (r)
            if indvEvacTime == -1 % no decision; ALL P_hits seen until end of trial (P_hit = 1 or 0)
                assert(looPhits(end,j,n) == 0 || looPhits(end,j,n) == 1);
                h1 = histcounts(looPhits(1:find(looPhits(:,j,n)==looPhits(end,j,n),1,'first'),j,n),bins);
            else
                assert(looPhits(indvEvacTime,j,n)==ePhits(i,j,n));
                h1 = histcounts(looPhits(1:indvEvacTime,j,n),bins);
            end
            H(i,:,n) = H(i,:,n)+h1;
        end
    end
end

for n = 1:length(z)
    for i = 1:50 % iterate over all individials
        J(i,:,n) = J(i,:,n) + histcounts(ePhits(i,:,n),bins);
    end
end

for i = 1:50 % iterate through each individual
    for j = 1:length(z) % iterate through each trial
        indvEvacTime = evacTimes(i,j); % this is the time the individual (j) evacuated in trial (r)
        if indvEvacTime == -1 % no decision; ALL P_hits seen until end of trial (P_hit = 1 or 0)
            assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
            h1 = histcounts(rP_hits(1:find(rP_hits(:,j)==rP_hits(end,j),1,'first'),j),bins);
        else
            assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
            h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
        end
        Htrue(i,:) = Htrue(i,:)+h1;
    end
end

for i = 1:50 % iterate over all individials
    Jtrue(i,:) = Jtrue(i,:) + histcounts(evacPhits(i,:),bins);
end


H = sum(H,1);
J = sum(J,1);
H = H(:,2:end,:);
J = J(:,2:end,:);
Htrue = sum(Htrue,1);
Htrue = Htrue(2:end);
Jtrue = sum(Jtrue,1);
Jtrue = Jtrue(2:end);
P_hit_range = P_hit_range(2:end);


theta_p = [0.1; 0.1];
powerfun = @(theta_p)powerf(Htrue',Jtrue',theta_p,P_hit_range);
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[theta_power_all,~] = fminunc(powerfun,theta_p,options);
power_model = @(Ph,theta_p) theta_p(1)*Ph.^theta_p(2);

theta_cv=zeros(2,length(z));
for i = 1:length(z)
    theta = [0.1; 0.1];
    powerfun = @(theta_p)powerf(H(:,:,i),J(:,:,i),theta_p,P_hit_range); % cannot include P_hit = 0

    options = optimoptions(@fminunc,'MaxFunEvals',10000);
    [thetas,~] = fminunc(powerfun,theta,options);
    theta_cv(:,i) = thetas;
end
%%

endTimes = zeros(size(rP_hits,2),1);
for i = 1:size(rP_hits,2)
    endTimes(i) = find(rP_hits(:,i)==rP_hits(end,i),1,'first');
end

P_LOO_power = zeros(60,length(z));
Probs_LOO_power = cell(length(z),1);
mProbs_LOO_power = cell(length(z),1);
stdevs_LOO_power = cell(length(z),1);
T_LOO_power = cell(length(z),1);
for i = 1:length(z)
    [Ptest,Ttest,PPtest] = mastereq(power_model(rP_hits(:,i),theta_cv(:,i)),endTimes(i));
    mProbs_LOO_power{i} = Ptest;
    T_LOO_power{i} = Ttest;
    temp = interp1(Ttest,Ptest,1:1:endTimes(i));
    temp(endTimes(i)+1:60)=temp(end);    
    P_LOO_power(:,i) = temp;
    Probs_LOO_power{i} = PPtest;
    stdevs_LOO_power{i} = zeros(size(Ttest));
    for j = 1:length(Ttest)
        stdevs_LOO_power{i}(j) = sqrt(((0:50)-Ptest(j)).^2*(PPtest(j,:)'));
    end
end


%%
rss_cv = zeros(16,1);
% rss_cv_model = zeros(16,1);
for i = 1:length(z)
    for j = 1:endTimes(i)
        rss_cv(i) = rss_cv(i) + (P_LOO_power(j,i)-evac(j,i))^2;
%         rss_cv_model(i)=rss_cv_model(i)+(P_LOO_power(j,i)-P_power(j,i))^2;
    end
    rss_cv(i)=rss_cv(i)/endTimes(i);
%     rss_cv_model(i)=rss_cv_model(i)/60;
end
mse_cv = mean(rss_cv);
% mse_cv_model = mean(rss_cv_model);
rmse_cv = sqrt(mse_cv);

% save workspace
clear bins temp i j
save crossval_Ind50