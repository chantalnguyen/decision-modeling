%       evac    no evac
% hit     -6        -10
% miss    -2          0

hitEvac = -6;
hitNoEvac = -10;
missEvac = -2;
missNoEvac = 0;

% import data 
load('evacuate_data.mat');

gameinfo([30 74 121],end) = 1;

uniquePhits = unique(gameinfo,'rows');
endTimes = zeros(size(uniquePhits,1),1);
for i = 1:size(uniquePhits,1)
    endTimes(i) = find(uniquePhits(i,:)==uniquePhits(i,end),1,'first');
end

rPhits = 10.^floor(log10(abs(uniquePhits)));
rPhits = floor(uniquePhits./rPhits).*rPhits;
rPhits(isnan(rPhits)) = 0;
rPhits = round(rPhits,1);
%%
% define bounds for parameter space
alpha_lb = 0.5;
alpha_ub = 1.0;
beta_lb = 4;
beta_ub = 10;

% number of grid elements in one dimension
ngrid = 10;

% define grid of parameters
[grid_a,grid_b] = meshgrid(linspace(alpha_lb,alpha_ub,ngrid),linspace(beta_lb,beta_ub,ngrid));

power_model = @(Ph,theta) theta(1)*Ph.^theta(2);

endProbs = zeros(ngrid*ngrid,size(uniquePhits,1));

for g = 1:size(uniquePhits,1)
    count = 1;
    for j = 1:ngrid
        for i = 1:ngrid
            [~,~,endP] = indivmastereq(power_model(rPhits(g,1:endTimes(g)),[grid_a(i,j),grid_b(i,j)]),endTimes(g),1);
            endProbs(count,g) = endP;
            count = count+1;
        end
    end
end

scores = zeros(size(endProbs));
for g = 1:size(uniquePhits,1)
    if uniquePhits(g,end) == 1
        scores(:,g) = endProbs(:,g)*hitEvac + (1-endProbs(:,g))*hitNoEvac;
    else
        scores(:,g) = endProbs(:,g)*missEvac + (1-endProbs(:,g))*missNoEvac;
    end
end
%%
endscores = sum(scores,2);
endscores = endscores + 10*size(uniquePhits,1);
[maxscore,I]=max(endscores);
[r,c]=ind2sub(size(grid_a),I);
max_a = grid_a(r,c);
max_b = grid_b(r,c);
