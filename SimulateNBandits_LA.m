function data = SimulateNBandits_LA(task, params)
%% LAZY ARROWS 
% pick side most arrows are pointing to
%
% data - matrix with Ntrials rows
% [t cb iter stim rew_incorr rew_corr b s corr r prob]
% t: trial number. 1->Ntrials
% cb: correct bandit- 1->Nbandits
% iter: iteration number in the block. 1->?
% stim: Stimulus [0:Ndirections-1]^Nbandits
% rew_incorr: reward on trial if incorrect response (0 or 1)
% rew_corr: reward on trial if correct response (0 or 1)
% b: chosen bandit (unobserved by experimenter)
% s: chosen side (observed)
% corr: correct? (unobserved by participant) 
% r: reward? (observed) 
% prob: probability of choosing each bandit based on overlap

%% set task params
Ntrials = task.Ntrials; %number of trials in the session
stimData = task.stimData; %Stimulus Data [t cb iter stim rew_incorr rew_corr]

%% set model params
%default 
epsilon = 0;

if nargin > 1
    if ~isempty(params)
        epsilon = params(1); %noise
    end
end

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

    %introduce noise after selecting bandit
    if rand<epsilon
        if Ndirections == 2 %only L/R
            s = 1-s; %switch to opposite side
        else
            opts = setdiff(1:Ndirections-1,s);
            s = randsample(opts); %randomly sample from other direction options
        end
    end
    
    % correct side?
    corr = s==stim(cb);
    % probabilistic reward
    r = rew(1+corr);
    
    model_data = [model_data;[b s corr r prob]];
         
end

data = [stimData, model_data];
end
