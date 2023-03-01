function data = get_stim(task, visualize)
% generates stimuli sequence for one session
% option to use weighting on categories of stimuli combinations:
% all-same, all-unique, and some-different

if nargin == 1
    visualize = 0;
end

if isfield(task, 'weight')
    weight = task.weight;
else
    weight = 0;
end

if weight
    sud_ratio = task.sud_ratio;
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

%find indices & counts of types of trials
all_same = find(counts == Nbandits);
all_unique = find(counts == 1);
some_diff = find(counts ~= Nbandits & counts ~= 1);
Nsame = length(all_same);
Nunique = length(all_unique);
Ndiff = length(some_diff);

%get probabilities of each combination
if ~weight
    %all equal weighting, using natural statistics based on task parameters
    pcombos = ones(Ncombos, 1)/length(combos);
else
    %use input ratio of all-same to all-unique to some-diff
    s = sud_ratio(1); u = sud_ratio(2); d = sud_ratio(3);

    if Nunique || u == 0
        punique = (u/(sum(sud_ratio)))/Nunique;
        psame = (s/(sum(sud_ratio)))/Nsame;
        pdiff = (d/(sum(sud_ratio)))/Ndiff;
    else
        %warning('U set to non-zero but no available unique sets')
        punique = 0;

        %reweight same and some different
        psame = (s/(s+d))/Nsame;
        pdiff = (d/(s+d))/Ndiff;
    end

    pcombos = nan(Ncombos, 1);
    pcombos(all_same) = psame;
    pcombos(all_unique) = punique;
    pcombos(some_diff) = pdiff;
end


%% generate stimuli sequences
stim_inds = NaN(Ntrials, 1);

data = [];
cb = 1; %always start correct bandit as bandit 1
iter = 1;

% loop through trials
for t = 1:Ntrials
    % choose stim combo for trial
    stim_inds(t) = randsample(1:Ncombos, 1, true, pcombos); 
    stim = combos{stim_inds(t)};
    
    % reward if incorrect/correct (0 or 1)
    rew_incorr = rand<prew(1);
    rew_corr = rand<prew(2);
    
    data = [data;[t cb iter stim rew_incorr rew_corr]];
    
    if iter>10 & rand<pswitch %switch correct bandit
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
end

if visualize %visualize when testing only, do not recommend otherwise
    % figures will show distribution of overlap & categories compared to 
    % natural weighting

    dt = data(:, 4:3+Nbandits); %stimuli presentations in session
    dt = num2cell(dt, 2);
    c = cellfun(fun, dt); %max num overlap
    
    nat_overlap_count = hist(counts, 1:Nbandits); %natural weight overlap counts
    session_overlap_count = hist(c, 1:Nbandits); %session overlap counts
    
    figure
    bar([nat_overlap_count'./Ncombos, session_overlap_count'./Ntrials])
    ylim([0, 1])
    ylabel('Proportion')
    xlabel('Max Number of Arrows in Same Direction')
    legend({'Natural Weights', 'Actual Stimuli'})
    title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections))})

    %get category counts
    all_same_dt = find(c == Nbandits);
    all_unique_dt = find(c == 1);
    some_diff_dt = find(c ~= Nbandits & c ~= 1);
    Nsame_dt = length(all_same_dt);
    Nunique_dt = length(all_unique_dt);
    Ndiff_dt = length(some_diff_dt);

    nat_cat_counts = [Nsame, Nunique, Ndiff];
    session_cat_counts = [Nsame_dt, Nunique_dt, Ndiff_dt];
    
    figure
    bar([nat_cat_counts'./Ncombos, session_cat_counts'./Ntrials])
    ylim([0, 1])
    ylabel('Proportion')
    xlabel('Category')
    xticklabels({'All Same', 'All Unique', 'Some Different'})
    legend({'Natural Weights', 'Actual Stimuli'})
    title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections)), strcat('ratio = ', num2str(sud_ratio))})
end


end