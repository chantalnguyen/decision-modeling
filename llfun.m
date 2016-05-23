% Calculates (negative) log-likelihood function to be minimized, according
% to beta distribution of successes J and failures H-J

function ll = llfun(H,J,theta,P_hit_range,space)
if nargin < 5
    space = 50;
end
ll = 0;
if space == 50
    for i = 1:length(P_hit_range)
        ll = ll + (H(i) - J(i))*log(1 - theta(1)*P_hit_range(i).^theta(2)) + J(i)*log(theta(1)*...
            P_hit_range(i).^theta(2));
    end
else
    for i = 1:length(P_hit_range)
        ll = ll + (H(i) - J(i))*log(1 - theta(1)*sign(P_hit_range(i)).*abs(P_hit_range(i)^(theta(2)/space)))...
            + J(i)*log(theta(1)*sign(P_hit_range(i)).*abs(P_hit_range(i)^(theta(2)/space)));
    end
end
ll = -1*ll;
end
