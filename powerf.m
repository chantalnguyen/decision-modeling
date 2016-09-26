
% Calculates (negative) log-likelihood function to be minimized, according
% to beta distribution of successes J and failures H-J

function power_ll = powerf(H,J,theta,P_hit_range)
power_ll = 0;
for i = 1:length(P_hit_range)
    power_ll = power_ll + (H(i) - J(i))*log(1 - theta(1)*P_hit_range(i).^theta(2)) + J(i)*log(theta(1)*...
        P_hit_range(i).^theta(2));
end
power_ll = -1*power_ll;
end