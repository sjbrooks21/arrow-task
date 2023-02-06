function data = SimulateNBandits_WSLS_side(task,params)
%% 
% Pick side and stick with side until no reward, then switch
% data - Ntrials * (7+2*Nbandits) matrix
%[t cb iter stim s cor r]];
% 1. trial number. 1->Ntrials
% 2. correct bandit- 1->Nbandits
% 3. iteration number in the block. 1->?
% 3+[1:nbandits]. Stimulus [0,1]^Nbandits 
% %%3+nbandits + 1 - chosen bandit (unobserved) (not choosing bandit tho)
% 3+nbandits + 1 - chosen side (observed)
% 3+nbandits + 2 - correct? (unobseverd) 
% 3+nbandits + 3 - reward? (obseverd) 

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
s = randi([0,1]); %initialize side choice, equal probability left vs right

for t = 1:Ntrials
    stim = randi([0,1],[1,Nbandits]);
    
    %introduce noise
    if rand<epsilon
        s = 1-s;
    end
    
    %bandit
    %b_prob = (stim == s)/(sum(stim == s));
   
    %b = find(mnrnd(1,b_prob));
    
    % correct side?
    cor = s==stim(cb);
    % probabilistic reward
    r = (rand<prew(1+cor));
    
    data = [data;[t cb iter stim s cor r]];
    
    %update next side based on reward
    if ~r %if no reward, switch side
        s = 1 - s;
    end
    
    if iter>10 && rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
        
    
end

