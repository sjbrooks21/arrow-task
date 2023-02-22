function data = SimulateNBandits_stick_track_correct(task,params)
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
nOut = (8+2*Nbandits);
cb = 1;
t_iter = 1;
c_iter = 0;
prior = ones(1,Nbandits)/Nbandits;
data = NaN(Ntrials, nOut);
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
            opts = setdiff((1:Ndirections) - 1,s);
            s = randsample(opts); %randomly sample from other direction options
        end
    end
    % correct side?
    cor = s==stim(cb);
    c_iter = c_iter + cor; %update correct count
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

