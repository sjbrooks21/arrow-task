function data = get_stim(task, weight, ref_dur, sd_ratio, visualize)

if nargin < 5
    visualize = 0;
end
if nargin < 4
    sd_ratio = 1; %equally likely to pick one of same vs one of different
end
if nargin < 3
    ref_dur = 1;
end
if nargin == 1
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

%find where all arrows in same direction
all_same = find(counts == Nbandits);
Nsame = length(all_same);
Ndiff = Ncombos-Nsame;

if weight
    pcombos = 1./counts; %weight inversely to number of arrows pointing in same direction
    psame = 1/Nbandits;
else
    %may change sd_ratio
    %determine ratio of psame to pdiff
    %sd_ratio = 0.5; %half as likely to pick one of same vs one of different
    %sd_ratio = 1; %equally likely to pick one of same vs one of different
    pdiff = 1/(Ndiff + sd_ratio*Nsame);
    psame = pdiff*sd_ratio;
    pcombos = ones(Ncombos, 1)*pdiff;
    pcombos(all_same) = psame;
end


stim_inds = NaN(Ntrials, 1);
num_same = 0;
max_same = floor(0.05*Ntrials); %max 5% of trials all same
same_ref = 0;

% %dictate how often you allow all arrows in same direction
% psame = 0.05;
% pdiff = 0.95;
% 
% %assign individual probabilities
% pcombos = ones(Ncombos, 1)*pdiff/Ndiff;
% pcombos(all_same) = psame/Nsame;

data = [];
cb = 1;
iter = 1;


for t = 1:Ntrials
    %stim = randi([0,Ndirections-1],[1,Nbandits]);
    stim_inds(t) = randsample(1:Ncombos, 1, true, pcombos);
    stim = combos{stim_inds(t)};
    
    %if picked all same, start refractory
    if ismember(stim_inds(t), all_same)
        num_same = num_same + 1;
        same_ref = ref_dur;
        pcombos(all_same) = 0;
        
    %if refractory, reduce
    elseif same_ref
        same_ref = same_ref - 1;
        
        %if no longer refractory & under max, reset probs to psame
        if ~same_ref && num_same < max_same
            pcombos(all_same) = psame;  
        end
    end
    
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