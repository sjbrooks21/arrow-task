function data = get_stim3(task, visualize)

if nargin == 1
    visualize = 0;
end

%get task parameters
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
prew = task.prew;

%get all possible stimuli presentations
combos = num2cell(permn(0:Ndirections - 1, Nbandits), 2);
Ncombos = length(combos);

%get max number of same direction
fun = @(x) max(hist(x, 0:Ndirections-1));
counts = cellfun(fun, combos);

%find where count of types of trials
all_same = find(counts == Nbandits);
all_unique = find(counts == 1);
some_diff = find(counts ~= Nbandits & counts ~= 1);
Nsame = length(all_same);
Nunique = length(all_unique);
Ndiff = length(some_diff);

%all equal weighting, using natural statistics based on task parameters
pcombos = ones(Ncombos, 1)/length(combos);

stim_inds = NaN(Ntrials, 1);


data = [];
cb = 1;
iter = 1;


for t = 1:Ntrials
    %stim = randi([0,Ndirections-1],[1,Nbandits]);
    stim_inds(t) = randsample(1:Ncombos, 1, true, pcombos);
    stim = combos{stim_inds(t)};
    
    rew_incorr = rand<prew(1);
    rew_corr = rand<prew(2);
    
    data = [data;[t cb iter stim rew_incorr rew_corr]]; %trial, correct bandit, iteration, stimuli, reward
    
    if iter>10 & rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
end

if visualize
    dt = data(:, 4:3+Nbandits);
    dt = num2cell(dt, 2);
    c = cellfun(fun, dt);
    
    c1 = hist(counts, 1:Nbandits);
    c2 = hist(c, 1:Nbandits);
    
    figure
    bar([c1'./Ncombos, c2'./Ntrials])
    ylim([0, 1])
    ylabel('Proportion')
    xlabel('Max Number of Arrows in Same Direction')
    legend({'Combos', 'Actual Stimuli'})
    title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections))})

%     s = ismember(stim_inds, all_same);
%     figure;
%     h = stem(1:Ntrials, s);
%     set(h, 'Marker', 'none')
%     title({'Trials with all same direction', 'semi-controlled', strcat('ref =', num2str(ref_dur), ' weighted = ', num2str(weight))})
%     xlabel('trial')
%     ylabel('is all same')
end


end