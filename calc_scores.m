% Calculate scores
% Input arguments:
%   evac_data: an array containing the columns of peoeva corresponding to
%       the desired trials
%   trial_numbers: a vector of the desired trial numbers
%   gameinfo: the array containing the P_hit trajectories of all trials
% Output variables:
%   scoremat: an array containing the scores of each participant for each
%       of the desired trials (scoremat is the same size as evac_data)

function scoremat = calc_scores(evac_data,trial_numbers,gameinfo)
scoremat = evac_data;
    for i = 1:length(trial_numbers)
        if gameinfo(trial_numbers(i),end) == 1 || trial_numbers(i) == 30 || trial_numbers(i) == 74 || trial_numbers(i) == 121 % hit
            scoremat(:,i) = scoremat(:,i)*4;
        else % miss
            scoremat(:,i) = scoremat(:,i)*-2 + 10;
        end
    end
end