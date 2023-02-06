function data = SimulateNBandits_LA(task)
%% 
% LAZY ARROWS 
% pick side most arrows are pointing to
% data - Ntrials * (7+Nbandits) matrix
%[t cb iter stim b s cor r]];
% 1. trial number. 1->Ntrials
% 2. correct bandit- 1->Nbandits
% 3. iteration number in the block. 1->?
% 3+[1:nbandits]. Stimulus [0,Ndirections]^Nbandits 
% 3+nbandits + 1 - chosen bandit (unobserved) (not choosing bandit tho)
% 3+nbandits + 2 - chosen side (observed)
% 3+nbandits + 3 - correct? (unobseverd) 
% 3+nbandits + 4 - reward? (obseverd) 

%% set task params
prew = task.prew;
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
stimData = task.stimData; %[t cb iter stim rew_incorr rew_corr]

%%

model_data = [];
for t = 1:Ntrials
    cb = stimData(t, 2);
    stim = stimData(t, 4:end-2);
    rew = stimData(t, end-1:end);
    
    % side choice
    [M,~,C] = mode(stim);
    if length(C) > 1 %there is a tie
        s = randsample(C{1}); %randomly choose one of sides that tied
    else
        s = M;
    end
    
    %bandit 
    %(not really relevant here since not choosing side based on bandit, but
    % selecting one of bandits anyway) 
    prob = (stim == s)/(sum(stim == s)); %equal nonzero prob assigned to bandits pointing in chosen direction
   
    b = find(mnrnd(1,prob));
    
    % correct side?
    cor = s==stim(cb);
    % probabilistic reward
    r = rew(1+cor);
    
    model_data = [model_data;[b s cor r prob]];
         
end

data = [stimData, model_data];
end
