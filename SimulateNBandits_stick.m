function data = SimulateNBandits_stick(task,params)
%% set task params

prew=task.prew;
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
stimData = task.stimData;

%% set model params
beta = params(1);
epsilon = params(2);
stick = params(3);

%%

prior = ones(1,Nbandits)/Nbandits;
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
    model_data = [model_data;[b s cor r prob]];
        
    
end

data = [stimData, model_data];
end
