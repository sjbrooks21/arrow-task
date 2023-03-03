function data = SimulateNBandits_RL(task,params)
%% Skylar Brooks 03/03/2023
% rotation project in CCN lab
% extended from previous code by Anne Collins- CCN Lab
%
%% Reinforcement Learning
% Choose bandits based on Q-values, update Q-values after reward outcome
% Has additional stick parameter which increases likelihood of staying with
% bandit last chosen
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
alpha = params{1}; %learning rate array, [alpha_minus, alpha_plus]
stick = params{2}; %choice stickiness

%defaults
beta = 20; %stochasticity
epsilon = 0; %noise

%override defaults if included in parameter input
if length(params) > 2
    beta = params{3};    

    if length(params) > 3
        epsilon = params{4};
    end
end
    
%%

Q = ones(1,Nbandits)/Nbandits; %initialize Qs to have equal values
model_data = [];

% loop over trials
for t = 1:Ntrials
    cb = stimData(t, 2);
    stim = stimData(t, 4:end-2);
    rew = stimData(t, end-1:end);
    
    % bandit choice
    W = Q;
    if t>1
        W(b) = W(b)+stick;
    end
    
    % compute the softmax probability
    prob = exp(beta*W);
    prob = prob/sum(prob);
    
    % make a choice
    b = find(mnrnd(1,prob)); 
    
    % side choice
    s = stim(b);
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
    
    % update arrow value - RL delta rule
    Q(b) = Q(b) + alpha(r+1) * (r-Q(b)); %if r, Q increases, else Q decreases
    %alpha_minus used if r = 0, alpha_plus used if r = 1
    
    % counterfactual updating
    others = setdiff(1:Nbandits,b); %find indices of arrows not chosen
    Q(others) = Q(others) + alpha(r+1)* ((1-r)-Q(others));
    
    % storing the data
    model_data = [model_data;[b s corr r prob]];
        
    
end

data = [stimData, model_data];
end
