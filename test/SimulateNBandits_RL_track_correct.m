function data = SimulateNBandits_RL_track_correct(task,params)
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

%% set model params
beta = 20;%params(1); %this is stochasticity (randomness) param
epsilon = 0;%params(2); %this is noise
alpha = params(1); %learning rate
stick = params(2); %choice stickiness

%%
nOut = (8+2*Nbandits);

cb = 1;
t_iter = 1;
c_iter = 0;
Q = ones(1,Nbandits)/Nbandits; %initialize Qs to have equal values
data = NaN(Ntrials, nOut);

% loop over trials
for t = 1:Ntrials
    stim = randi([0,1],[1,Nbandits]);
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

    %introduce noise
    if rand<epsilon
        if Ndirections == 2 %only L/R
            s = 1-s; %switch to opposite side
        else
            opts = setdiff((1:Ndirections)-1,s);
            s = randsample(opts); %randomly sample from other direction options
        end
    end
    
    % correct side?
    cor = s==stim(cb);
    c_iter = c_iter + cor; %update correct count
    % probabilistic reward
    r = (rand<prew(1+cor)); %if cor, choose prew(2) which is prew if correct
    
    % update arrow value - RL delta rule
    Q(b) = Q(b) +alpha * (r-Q(b)); %if r, Q increases, else Q decreases
    
    % counterfactual updating
    others = setdiff(1:Nbandits,b); %find indices of arrows not chosen
    Q(others) = Q(others) + alpha* ((1-r)-Q(others));
    
    % storing the data
    data(t, :) = [t cb t_iter c_iter stim b s cor r Q];
    
    % change in correct arrow after at least 10 trials
    if c_iter>=10 & rand<pswitch
        t_iter = 1;
        c_iter = 0;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        t_iter = t_iter + 1;
    end
        
    
end

