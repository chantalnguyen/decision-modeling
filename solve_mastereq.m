%% solve master equation
power_model = @(Ph,theta_i) theta_i(1)*Ph.^theta_i(2);

% calculate endtimes of each trial
endTimes = zeros(size(rP_hits,2),1);
for i = 1:size(rP_hits,2)
    endTimes(i) = find(rP_hits(:,i)==rP_hits(end,i),1,'first');
end
theta = [1, 5.899];
P_statopt = zeros(60,length(z)); % mean cumulative evacuations for each trial
Probs_statopt = cell(length(z),1); % full probability distributions for each trial
mProbs_statopt = cell(length(z),1); % (raw) mean cumulative evacuations for each trial
stdevs_statopt = cell(length(z),1); % standard deviation of mean cumulative evacuations
T_statopt = cell(length(z),1); % ODE time steps for each trial
for i = 1:length(z)
    [Ptest,Ttest,PPtest] = mastereq(power_model(rP_hits(:,i),theta),endTimes(i));
    mProbs_statopt{i} = Ptest;
    T_statopt{i} = Ttest;
    temp = interp1(Ttest,Ptest,1:1:endTimes(i));
    temp(endTimes(i)+1:60)=temp(end);
    P_statopt(:,i) = temp;
    Probs_statopt{i} = PPtest;
    stdevs_statopt{i} = zeros(size(Ttest));
    for j = 1:length(Ttest)
        stdevs_statopt{i}(j) = sqrt(((0:50)-Ptest(j)).^2*(PPtest(j,:)'));
    end
end

% calculate root mean square error
rss_statopt = zeros(length(z),1);
for i=1:length(z)
    for j = 1:endTimes(i)
        rss_statopt(i)=rss_statopt(i)+(P_statopt(j,i)-evac(j,i))^2;
    end
    rss_statopt(i)=rss_statopt(i)/endTimes(i);
end
mse_statopt = mean(rss_statopt);
rmse_statopt = sqrt(mse_statopt);

% save workspace
clear temp i j
