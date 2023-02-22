function data = SimulateNBandits_WSLS_track_correct(task, params)
%% 
% Pick arrow and stick with arrow until no reward, then switch to arrow
% pointing in another direction
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
Ndirections = task.Ndirections;

%% set model params
epsilon = params(1);

%%
nOut = (8+2*Nbandits);
cb = 1;
t_iter = 1;
c_iter = 0;
data = NaN(Ntrials, nOut);
prob = ones(1,Nbandits)/Nbandits;

for t = 1:Ntrials
    stim = randi([0,Ndirections-1],[1,Nbandits]);
    
    %bandit
    b = randsample(1:Nbandits, 1, true, prob);
    
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
    r = (rand<prew(1+cor));
        
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
    
    % storing the data
    data(t, :) = [t cb t_iter c_iter stim b s cor r prob];
    
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

