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

prew=task.prew; %probability of reward
pswitch = task.pswitch; %probability of switching
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;

%% set model params
epsilon = 0;%params(2); %this is noise
alpha = params(1); %learning rate
stick = params(2); %choice stickiness

if length(params) > 2
    beta = params(3);
else
    beta = 20;
end
    
%%
nOut = (7+2*Nbandits);

cb = 1;
iter = 1;
Q = ones(1,Nbandits)/Nbandits; %initialize Qs to have equal values
data = NaN(Ntrials, nOut);

% loop over trials
for t = 1:Ntrials
    stim = randi([0, Ndirections-1],[1,Nbandits]);
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
    r = (rand<prew(1+cor)); %if cor, choose prew(2) which is prew if correct
    
    % update arrow value - RL delta rule
    Q(b) = Q(b) +alpha * (r-Q(b)); %if r, Q increases, else Q decreases
    
    % counterfactual updating
    others = setdiff(1:Nbandits,b); %find indices of arrows not chosen
    Q(others) = Q(others) + alpha* ((1-r)-Q(others));
    
    % storing the data
    data(t, :) = [t cb iter stim b s cor r Q];
    
    % change in correct arrow after at least 10 trials
    if iter>10 & rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
        
    
end

