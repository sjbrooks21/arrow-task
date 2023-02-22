function data = SimulateNBandits_LA_track_correct(task)
%% 
% LAZY ARROWS 
% pick side most arrows are pointing to
% data - Ntrials * (7+Nbandits) matrix
%[t cb iter stim b s cor r]];
% 1. trial number. 1->Ntrials
% 2. correct bandit- 1->Nbandits
% 3. iteration number in the block. 1->?
% 3+[1:nbandits]. Stimulus [0,Ndirections]^Nbandits 
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

epsilon = 0;
%%
nOut = (8+2*Nbandits);
cb = 1;
t_iter = 1;
c_iter = 0;
data = NaN(Ntrials, nOut);
for t = 1:Ntrials
    stim = randi([0, Ndirections-1],[1,Nbandits]);
    
    % side choice
    [M,~,C] = mode(stim);
    if length(C) > 1 %there is a tie
        s = randsample(C{1}); %randomly choose one of sides that tied
    else
        s = M;
    end
    
    %introduce noise
    if rand<epsilon
        if Ndirections == 2 %only L/R
            s = 1-s; %switch to opposite side
        else
            opts = setdiff((1:Ndirections)-1,s);
            s = randsample(opts); %randomly sample from other direction options
        end
    end

    %bandit 
    %(not really relevant here since not choosing side based on bandit, but
    % selecting one of bandits anyway) 
    prob = (stim == s)/(sum(stim == s)); %equal nonzero prob assigned to bandits pointing in chosen direction
   
    b = find(mnrnd(1,prob));
    
    % correct side?
    cor = s==stim(cb);
    c_iter = c_iter + cor; %update correct count
    % probabilistic reward
    r = (rand<prew(1+cor));
    
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

