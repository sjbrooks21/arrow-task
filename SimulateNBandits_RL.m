function data = SimulateNBandits_RL(task,params)
%% 
% data - Ntrials * (7+2*Nbandits) matrix
%[t cb iter stim b s cor r Q]];
% 1. trial number. 1->Ntrials
% 2. correct bandit- 1->Nbandits
% 3. iteration number in the block. 1->?
% 3+[1:nbandits]. Stimulus [0,1]^Nbandits 
% 3+nbandits + 1 - chosen bandit (unobserved)
% 3+nbandits + 2 - chosen side (observed)
% 3+nbandits + 3 - correct? (unobseverd) 
% 3+nbandits + 4 - reward? (obseverd) 
% 3+nbandits + 4 + [1:nbandits] - bandit Q-values
%% set task params

prew = task.prew; %probability of reward
pswitch = task.pswitch; %probability of switching
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
stimData = task.stimData; %[t cb iter stim rew_incorr rew_corr]

%% set model params
epsilon = 0; %noise
alpha = params{1}; %learning rate
stick = params{2}; %choice stickiness

if length(params) > 2
    beta = params{3};
else
    beta = 20;
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
    cor = s==stim(cb);
    % probabilistic reward
    r = rew(1+cor);
    
    % update arrow value - RL delta rule
    Q(b) = Q(b) +alpha(r+1) * (r-Q(b)); %if r, Q increases, else Q decreases
    
    % counterfactual updating
    others = setdiff(1:Nbandits,b); %find indices of arrows not chosen
    Q(others) = Q(others) + alpha(r+1)* ((1-r)-Q(others));
    
    % storing the data
    model_data = [model_data;[b s cor r prob]];
        
    
end

data = [stimData, model_data];
end
