function X = forfit(data,task)

%X: stimulus, side, reward

% data - Ntrials * (7+2*Nbandits) matrix
%[t cb iter stim b s cor r Q]];
% 1. trial number. 1-Ntrials
% 2. correct bandit- 1-Nbandits
% 3. iteration number in the block. 1-?
% 3+[1:nbandits]. Stimulus [0,1]^Nbandits 
% 3+nbandits + 1 - chosen bandit (unobserved)
% 3+nbandits + 2 - chosen side (observed)
% 3+nbandits + 3 - correct? (unobseverd) 
% 3+nbandits + 4 - reward? (obseverd) 
% 3+nbandits + 4 + [1:nbandits] - bandit Q-values
Nbandits = task.Nbandits;

X = data(:,3+([1:Nbandits, Nbandits + [2 4]]));
end