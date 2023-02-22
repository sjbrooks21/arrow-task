function data = SimulateNBandits_stick(task,params)
%% set task params

prew=task.prew;
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;

%% set model params
beta = params(1);
epsilon = params(2);
stick = params(3);

%%

cb = 1;
iter = 1;
prior = ones(1,Nbandits)/Nbandits;
data = [];
for t = 1:Ntrials
    stim = randi([0,Ndirections-1],[1,Nbandits]);
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
    r = (rand<prew(1+cor));
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
    
    data = [data;[t cb iter stim b s cor r p]];
    
    if iter>10 & rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
        
    
end

