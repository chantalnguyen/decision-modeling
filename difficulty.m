load experiment
load peoeva
gameinfo([30 74 121],end)=1;
trialOutcome = gameinfo(:,end); % if the trial hit or miss
numEvac = sum(peoeva);
difficulty = numEvac;
difficulty(trialOutcome==1) = difficulty(trialOutcome==1)./50;
difficulty(trialOutcome==0) = (50-difficulty(trialOutcome==0))./50;
difficulty = difficulty';
difficulty(missing) = [];