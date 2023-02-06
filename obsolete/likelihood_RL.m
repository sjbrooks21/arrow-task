function data = likelihood_RL(params,X,Nbandits)
%% 
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
%% set task params

Nbandits = task.Nbandits;
stimuli = X(:,1:Nbandits);
choice = X(:,Nbandits+1);
reward = X(:,Nbandits+2);
%% set model params
beta = 20;%params(1);
epsilon = 0;%params(2);
alpha = params(1);
stick = params(2);

%%

Q = ones(1,Nbandits)/Nbandits;
data = [];
for t = 1:Ntrials
    stim = randi([0,1],[1,Nbandits]);
    % bandit choice
    W = Q;
    if t>1
        W(b) = W(b)+stick;
    end
    prob = exp(beta*W);prob = prob/sum(prob);
    b = find(mnrnd(1,prob));
    % side choice
    s = stim(b);
    if rand<epsilon
        s = 1-s;
    end
    % correct side?
    cor = s==stim(cb);
    % probabilistic reward
    r = (rand<prew(1+cor));
    % compute likelihood
    Q(b) = Q(b) +alpha * (r-Q(b));
    others = setdiff(1:Nbandits,b);
    Q(others) = Q(others) + alpha* ((1-r)-Q(others));
    
    data = [data;[t cb iter stim b s cor r Q]];
    
    if iter>10 & rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
        
    
end

