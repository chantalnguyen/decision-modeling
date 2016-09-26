% Perform Bayesian model selection for decision making as a function of
% disaster probability (with beta likelihood function)

% desired trial #s
z = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
z = z(randperm(length(z)));

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

% this is our data set
rP_hits = round_P_hits(z,:); % rounded P_hit trajectories for each trial
rP_hits = rP_hits'; 
evac = evapeocumu(z,:); % empirical cumulative evacuations for each trial
evac = evac';
evacTimes = evacuateTime(:,z); % times of evacuation decisions for each trial
evacPhits = evacuateProb(:,z); % P_hits of evacuation decisions for each trial

% round evacuateProb values down to nearest tenth
evacPhits2 = 10.^floor(log10(abs(evacPhits)));
evacPhits2 = floor(evacPhits./evacPhits2).*evacPhits2;
evacPhits2(isnan(evacPhits2))=0;
evacPhits = round(evacPhits2,1);
clear evacPhits2;

P_hit_range=0:0.1:1;

% Number of times participants at home observe P_hit values
H = zeros(length(P_hit_range),length(z)); 

% Number of times participants decide to evacuate at P_hit values
J = zeros(length(P_hit_range),length(z));

bins = -0.05:0.1:1.05;

for j = 1:length(z) % iterate through each trial
    for i = 1:50    % iterate through each individual
        indvEvacTime = evacTimes(i,j);
        if indvEvacTime == -1
            assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
            h1 = histcounts(rP_hits(1:find(gameinfo(z(j),:)==gameinfo(z(j),end),1,'first'),j),bins);
        else
            assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
            h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
        end
        H(:,j) = H(:,j) + h1';
    end
end

for i = 1:length(z)
    J(:,i) = J(:,i) + histcounts(evacPhits(:,i),bins)';
end


% initialize models
power_model = @(Ph,theta) theta(1)*Ph.^theta(2);

% define bounds for parameter space
alpha_lb = 0.3;
alpha_ub = 1.0;
beta_lb = 0;
beta_ub = 10;

% number of grid elements in one dimension
ngrid = 200;

% define grid of parameters
[grid_a,grid_b] = meshgrid(linspace(alpha_lb,alpha_ub,ngrid),linspace(beta_lb,beta_ub,ngrid));

% define uniform prior
prior = ones(ngrid,ngrid)/ngrid^2;

likelihoods = ones([size(prior) length(z)]);

for i = 1:length(z)
    for j = 1:ngrid
        for k = 1:ngrid
            for m = 2:length(P_hit_range)
                likelihoods(j,k,i) = likelihoods(j,k,i)*betapdf(power_model(P_hit_range(m),[grid_a(j,k),grid_b(j,k)]),J(m,i)+1,H(m,i)-J(m,i)+1);
            end
        end
    end
end

% perform Bayesian updates
posteriors = zeros(size(likelihoods));

posteriors(:,:,1) = prior.*likelihoods(:,:,1);
evidence = sum(sum(posteriors(:,:,1)));
posteriors(:,:,1) = posteriors(:,:,1)/evidence;

for i = 2:length(z)
    posteriors(:,:,i) = posteriors(:,:,i-1).*likelihoods(:,:,i);
    evidence = sum(sum(posteriors(:,:,i)));
    posteriors(:,:,i) = posteriors(:,:,i)/evidence;
end
%%
maxvals = zeros(length(z),1);
maxparams = zeros(length(z),2);
for i = 1:length(z)
    temp = posteriors(:,:,i);
    [M,I]=max(temp(:));
    [r,c]=ind2sub(size(grid_a),I);
    maxvals(i) = M;
    maxparams(i,:) = [grid_a(r,c) grid_b(r,c)];
end

%%
colors = flipud(jet(length(z)));
for i = 1:length(z)
    plot((0:0.01:1),maxparams(i,1)*(0:0.01:1).^(maxparams(i,2)),'color',colors(i,:))
    hold on
    
end
plot((0:0.01:1),maxparams(end,1)*(0:0.01:1).^(maxparams(end,2)),'color','b','linewidth',4)
% plot((0:0.01:1),0.8824*(0:0.01:1).^(5.4053),'color','g','linewidth',4,'linestyle',':')

%%
fig = figure;
surfc(grid_a,grid_b,prior,'edgecolor','none');
title ('posterior after 0 trials')
ylabel('b')
xlabel('a')
zlabel('probability')
zlim([0 0.01])
f = getframe(fig);
[im,map]=rgb2ind(f.cdata,256);
imwrite(im,map,'figures/bayesPosterior.gif', 'Loopcount',inf);
for i = 1:length(z)
    surfc(grid_a,grid_b,posteriors(:,:,i),'edgecolor','none');
    title (['posterior after ' num2str(i) ' trials'])
    annotation('textbox', [0.1,0.1,0.1,0.1],'String', ['trial ' num2str(z(i))]);
    ylabel('b')
    xlabel('a')
    zlabel('probability')
    zlim([0 0.01])
    f = getframe(fig);
    [im,map]=rgb2ind(f.cdata,256);
    imwrite(im,map,'figures/bayesPosterior.gif','WriteMode','append','DelayTime',0.2);
    clf(fig)
end
    