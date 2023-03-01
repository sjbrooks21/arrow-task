function data = SimulateNBandits_WSLS(task, params)
%% Win-Stay-Lose-Shift (Smart Shift)
% Pick arrow and stick with arrow until no reward, then switch to arrow
% pointing in another direction
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
% prob: probability of choosing each bandit at start of trial

%% set task params
Ntrials = task.Ntrials; %number of trials in the session
Nbandits = task.Nbandits; %number of bandits presented
Ndirections = task.Ndirections; %number of possible directions to choose from
stimData = task.stimData; %Stimulus Data [t cb iter stim rew_incorr rew_corr]

%% set model params
epsilon = params(1); %noise

%%

model_data = [];
prob = ones(1,Nbandits)/Nbandits; %initialize bandits to have equal probability

for t = 1:Ntrials
    cb = stimData(t, 2);
    stim = stimData(t, 4:end-2);
    rew = stimData(t, end-1:end);
    
    %bandit
    b = randsample(1:Nbandits, 1, true, prob);
    
    % side choice
    s = stim(b);
    
    %introduce noise
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
    
    %update probabilities based on reward
    if r %if reward, chose same bandit next time
        prob = zeros(1, Nbandits);
        prob(b) = 1;
    elseif any(stim ~= stim(b)) %choose another bandit if at least one pointing in different direction
        prob = stim ~= stim(b);  
        prob = prob/(sum(prob));
    else %if all pointing in same direction, randomly choose another bandit
        prob = ones(1, Nbandits);
        prob(b)= 0;
        prob = prob/(sum(prob));
    end
   
end

data = [stimData, model_data];

end

