function data = SimulateNBandits_stick(task,params)
%% Bayes Model (Sticky)
% Choose bandits based on priors, update priors after reward outcome
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
prew = task.prew; %probability of reward | incorrect/correct (e.g. [0.1, 0.9])
pswitch = task.pswitch; %probability of switching to a different rule
Ntrials = task.Ntrials; %number of trials in the session
Nbandits = task.Nbandits; %number of bandits presented
Ndirections = task.Ndirections; %number of possible directions to choose from
stimData = task.stimData; %Stimulus Data [t cb iter stim rew_incorr rew_corr]

%% set model params
beta = params(1); %stochasticity
epsilon = params(2); %noise
stick = params(3); %stickiness

%%
prior = ones(1,Nbandits)/Nbandits; %initialize bandits to have equal probability
model_data = [];
for t = 1:Ntrials
    cb = stimData(t, 2);
    stim = stimData(t, 4:end-2);
    rew = stimData(t, end-1:end);
    
    % add stickiness
    W = log(prior);
    if t>1
        W(b) = W(b)+stick;
    end
    % bandit policy
    prob = exp(beta*W);prob = prob/sum(prob);
    % bandit choice
    b = find(mnrnd(1,prob));
    % side choice
    s = stim(b);
    if rand<epsilon %introduce noise
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
    % compute likelihood
    for i=1:Nbandits
        likelihood(i) = prew(1+(stim(i)==s));
    end
    if r==0
        likelihood = 1-likelihood;
    end
    % compute posterior
    posterior = likelihood.*prior;
    p = posterior/sum(posterior);
    prior = (1-pswitch)*p+pswitch*(1-p)/(sum(1-p));
    
    % storing the data
    model_data = [model_data;[b s corr r prob]];
        
    
end

data = [stimData, model_data];
end
