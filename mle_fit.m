% Determine parameters of power-law decision model with maximum likelihood
% estimation
% Inputs:
%   z: one-dimensional vector containing numbers of trials
%   gameinfo: array of P_hit values for each trial
%   evapeocumu: array of cumulative evacuations for each trial
%   evacuateTime: array of decision times for each trial
%   evacuateProb: array of decision P_hits for each trial
%
% Optional arguments:
%   space: shelter space (default = 50)
%
% Outputs:
%   theta: 2-element vector containing best-fit parameters

function [theta] = mle_fit(z,gameinfo,evapeocumu,evacuateTime,evacuateProb,space)
if nargin <6
    space = 50;
end

% round P_hits down to nearest tenth, since participants only see discrete
% values of P_hit
round_P_hits = 10.^floor(log10(abs(gameinfo)));
round_P_hits = floor(gameinfo./round_P_hits).*round_P_hits;
round_P_hits(isnan(round_P_hits))=0;
round_P_hits = round(round_P_hits,1);

evacuateTime(evacuateTime == 0) = 1; % evacuating at 0 = evacuating at 1

rP_hits = zeros(60,length(z)); % rounded P_hit trajectories for each trial
evac = rP_hits; % empirical cumulative evacuations for each trial
evacTimes = zeros(50,length(z)); % times of evacuation DECISIONS for each trial
evacPhits = evacTimes; % P_hits of evacuation DECISIONS for each trial

% this is our data set
for i = 1:length(z)
    rP_hits(:,i) = round_P_hits(z(i),:);
    evac(:,i) = evapeocumu(z(i),:); 
    evacTimes(:,i) = evacuateTime(:,z(i));
    evacPhits(:,i) = evacuateProb(:,z(i));
end

P_hit_range=0:0.1:1;

% Number of times participants at home observe P_hit values
% H = zeros(length(P_hit_range),length(z)); 
H = zeros(50,length(P_hit_range));

% Number of times participants evacuate at P_hit values
% J = zeros(length(P_hit_range),length(z)); 
J = zeros(50,length(P_hit_range));

bins = -0.05:0.1:1.05;

if space == 50
    for i = 1:50 % iterate through each individual
        for j = 1:length(z) % iterate through each trial
            indvEvacTime = evacTimes(i,j); % this is the time the individual (i) evacuated in trial (j)
            if indvEvacTime == -1 % no decision; ALL rP_hits seen until end of trial (P_hit = 1 or 0)
%                 assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
                h1 = histcounts(rP_hits(1:find(rP_hits(:,j)==rP_hits(end,j),1,'first'),j),bins);
            else % decision; count P_hits seen until time of decision
%                 assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
                h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
            end
            H(i,:) = H(i,:)+h1;
        end
    end
else
    for i = 1:50 % iterate through each individual
        for j = 1:length(z) % iterate through each trial
            indvEvacTime = evacTimes(i,j); % this is the time the individual (i) evacuated in trial (j)
            if indvEvacTime == -1 % no decision; ALL rP_hits seen until when shelter is full NOT end of trial (P_hit = 1 or 0) 
%                 assert(rP_hits(end,j) == 0 || rP_hits(end,j) == 1);
                if evac(end,j) < space % if shelter did not fill, count all rP_hits seen
                    h1 = histcounts(rP_hits(1:find(rP_hits(:,j)==rP_hits(end,j),1,'first'),j),bins); 
                else % if shelter fills up, count rP_hits seen until shelter is full
                     % honestly, this would cause a discrepancy if
                     % participants decide when the shelter is already
                     % filled
                    h1 = histcounts(rP_hits(1:find(evac(:,j)==space,1,'first'),j),bins); 
                end
            else 
%                 assert(rP_hits(indvEvacTime,j)==evacPhits(i,j));
                h1 = histcounts(rP_hits(1:indvEvacTime,j),bins);
            end
            H(i,:) = H(i,:)+h1;
        end
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
[theta,~] = fminunc(powerfun,theta_0,options);

end
