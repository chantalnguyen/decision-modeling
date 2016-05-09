function [mean_evac,T,P] = mastereq(q,tf,ti,P_init,N)
% Evaluation of master equation for P_hit-dependent decision model
% Inputs:
%   q: decision model q (typically a power-law function evaluated for
%       specified parameters and P_hit trajectory) 
%   tf: final time step
%
% Optional arguments:
%   ti: initial time step
%   P_init: initial probability vector
%   N: total number of states (typically the number of individuals + 1)
%
% Outputs:
%   mean_evac: mean number of evacuations at each time step
%   T: vector of time steps
%   P: full probability distribution of possible states at each time step

tic;

if nargin<3
    
    ti = 1;
    N = 51;

    % initial probability distribution with all individuals at home
    P_init = zeros(N,1);
    P_init(1) = 1;
end

T_range = [ti tf];

%% solve master equation
[T, P] = ode45(@(t,P) odefunc(t,P,N,q), T_range, P_init, []);

% calculate mean number evacuated
mean_evac = P*(0:(N-1))';
toc;

%% the master equation
function dP=odefunc(t,P,N,q)
% initialize generator matrix
A = zeros(N,N); 

for n = 0:N-1
    for i = 0:N-1
        if i < n % can only transition from states of fewer to greater (or equal) evacuated
            qq = q(floor(t)); % round time down to nearest time step, since P_hit changes in discrete increments
            A(n+1,i+1)=nchoosek(N-i,n-i)*qq^(n-i)*(1-qq)^(N-n); % calculate transition rate
        end
    end
end
for i = 1:N
    A(i,i) = -sum(A(:,i)); % diagonal elements are negative sums of each column
end
 
dP=A*P;
end


end