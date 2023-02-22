function data = SimulateNBandits_WSLS_dumb(task, params)
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
prew = task.prew;
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
stimData = task.stimData; %[t cb iter stim rew_incorr rew_corr]

%% set model params
epsilon = params(1);

%%

model_data = [];
prob = ones(1,Nbandits)/Nbandits;

for t = 1:Ntrials
    cb = stimData(t, 2);
    stim = stimData(t, 4:end-2);
    rew = stimData(t, end-1:end);
    
    %bandit
    b = randsample(1:Nbandits, 1, true, prob);
    
    % side choice
    s = stim(b);
    
    %introduce noise
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
    
    model_data = [model_data;[b s cor r prob]];
    
    %update probabilities based on reward
    if r %if reward, chose same bandit next time
        prob = zeros(1, Nbandits);
        prob(b) = 1;
    else %randomly choose another bandit (does not track which one pointed in same direction)
        prob = ones(1, Nbandits);
        prob(b)= 0;
        prob = prob/(sum(prob));
    end
   
end

data = [stimData, model_data];

end

