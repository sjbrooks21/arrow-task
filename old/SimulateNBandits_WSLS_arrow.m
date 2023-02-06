function data = SimulateNBandits_WSLS_arrow(task, params)
%% 
% Pick arrow and stick with arrow until no reward, then switch
% data - Ntrials * (7+Nbandits) matrix
%[t cb iter stim b s cor r]];
% 1. trial number. 1->Ntrials
% 2. correct bandit- 1->Nbandits
% 3. iteration number in the block. 1->?
% 3+[1:nbandits]. Stimulus [0,1]^Nbandits 
% 3+nbandits + 1 - chosen bandit (unobserved) (not choosing bandit tho)
% 3+nbandits + 2 - chosen side (observed)
% 3+nbandits + 3 - correct? (unobseverd) 
% 3+nbandits + 4 - reward? (obseverd) 

%% set task params
prew=task.prew;
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;

%% set model params
epsilon = params(1);

%%

cb = 1;
iter = 1;
data = [];
prob = ones(1,Nbandits)/Nbandits;

for t = 1:Ntrials
    stim = randi([0,1],[1,Nbandits]);
    
    %bandit
    b = find(mnrnd(1,prob));
    
    % side choice
    s = stim(b);
    
    %introduce noise
    if rand<epsilon
        s = 1-s;
    end
    
    % correct side?
    cor = s==stim(cb);
    % probabilistic reward
    r = (rand<prew(1+cor));
    
    data = [data;[t cb iter stim b s cor r]];
    
    %update probabilities based on reward
    if r %if reward, chose same bandit next time
        prob = zeros(1, Nbandits);
        prob(b) = 1;
    else %else switch bandit
        prob = ones(1, Nbandits);
        prob(b) = 0; 
        prob = prob/(sum(prob)); %equal prob across other options
    end
    
    if iter>10 && rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
        
    
end

