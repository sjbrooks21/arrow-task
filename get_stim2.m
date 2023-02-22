function data = get_stim2(task, weight, sud_ratio, visualize)

if nargin == 3
    visualize = 0;
elseif nargin == 2
    sud_ratio = [1, 1, 1]; %equally likely to pick from any one of groups
elseif nargin == 1
    weight = 0;
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

if ~weight
    %all equal weighting, using natural statistics based on task parameters
    pcombos = ones(Ncombos, 1)/length(combos);
else
    %use input ratio of psame to punique to pdiff
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
    
    data = [data;[t cb iter stim rew_incorr rew_corr]];
    
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
    
%     c1 = hist(counts, 1:Nbandits);
%     c2 = hist(c, 1:Nbandits);
%     
%     figure
%     bar([c1'./Ncombos, c2'./Ntrials])
%     ylim([0, 1])
%     ylabel('Proportion')
%     xlabel('Max Number of Arrows in Same Direction')
%     legend({'Combos', 'Actual Stimuli'})
%     title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections))})

    all_same_dt = find(c == Nbandits);
    all_unique_dt = find(c == 1);
    some_diff_dt = find(c ~= Nbandits & c ~= 1);
    Nsame_dt = length(all_same_dt);
    Nunique_dt = length(all_unique_dt);
    Ndiff_dt = length(some_diff_dt);

    c12 = [Nsame, Nunique, Ndiff];
    c22 = [Nsame_dt, Nunique_dt, Ndiff_dt];
    
    figure
    bar([c12'./Ncombos, c22'./Ntrials])
    ylim([0, 1])
    ylabel('Proportion')
    xlabel('Category')
    xticklabels({'All Same', 'All Unique', 'Some Different'})
    legend({'Combos', 'Actual Stimuli'})
    title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections)), strcat('ratio = ', num2str(sud_ratio))})


%     s = ismember(stim_inds, all_same);
%     figure;
%     h = stem(1:Ntrials, s);
%     set(h, 'Marker', 'none')
%     title({'Trials with all same direction', 'semi-controlled', strcat('ref =', num2str(ref_dur), ' weighted = ', num2str(weight))})
%     xlabel('trial')
%     ylabel('is all same')
end


end