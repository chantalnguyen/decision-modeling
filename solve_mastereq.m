%% solve master equation
power_model = @(Ph,theta_i) theta_i(1)*Ph.^theta_i(2);

% calculate endtimes of each trial
endTimes = zeros(size(rP_hits,2),1);
for i = 1:size(rP_hits,2)
    endTimes(i) = find(rP_hits(:,i)==rP_hits(end,i),1,'first');
end
% theta = [1, 5.899];

P_bayes = zeros(60,length(z)); % mean cumulative evacuations for each trial
Probs_bayes = cell(length(z),1); % full probability distributions for each trial
mProbs_bayes = cell(length(z),1); % (raw) mean cumulative evacuations for each trial
stdevs_bayes = cell(length(z),1); % standard deviation of mean cumulative evacuations
T_bayes = cell(length(z),1); % ODE time steps for each trial
for i = 1:length(z)
    ind = find(z(i)==goodtrials);
    [Ptest,Ttest,PPtest] = mastereq(power_model(rP_hits(:,i),maxparams(ind-1,:)),endTimes(i));
    mProbs_bayes{i} = Ptest;
    T_bayes{i} = Ttest;
    temp = interp1(Ttest,Ptest,1:1:endTimes(i));
    temp(endTimes(i)+1:60)=temp(end);
    P_bayes(:,i) = temp;
    Probs_bayes{i} = PPtest;
    stdevs_bayes{i} = zeros(size(Ttest));
    for j = 1:length(Ttest)
        stdevs_bayes{i}(j) = sqrt(((0:50)-Ptest(j)).^2*(PPtest(j,:)'));
    end
end

% calculate root mean square error
rss_bayes = zeros(length(z),1);
for i=1:length(z)
    for j = 1:endTimes(i)
        rss_bayes(i)=rss_bayes(i)+(P_bayes(j,i)-evac(j,i))^2;
    end
    rss_bayes(i)=rss_bayes(i)/endTimes(i);
end
mse_bayes = mean(rss_bayes);
rmse_bayes = sqrt(mse_bayes);

% save workspace
clear temp i j
