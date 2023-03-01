addpath('..')
%% set up task parameters
task = struct;

% probability of reward|incorrect/correct left-right choice
task.prew=[.1 .9];
% probability of correct arrow changing
task.pswitch = .05;
% number of trials
task.Ntrials = 500;
% number of arrows
task.Nbandits = 3;
% number of observable actions (possible directions)
task.Ndirections = 4;
% weight of stimuli categories [0/1]
task.weight = 1;
% weight assignments
task.sud_ratio = [1, 1, 1];

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s corr r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

%%
Nbandits = [3, 4, 5]; %Nbandits to test
Ndirections = 4; %Ndirections to test

for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
            
            get_stim(task, 1); %includes visualization, change to 0 to skip
        end
end
