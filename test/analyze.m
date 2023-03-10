clear all

%% set up task parameters
task = [];

% probability of reward|incorrect/correct left-right choice
task.prew=[.1 .9];
%task.prew=[.1 .7];
% probability of correct arrow changing
task.pswitch = .05;
% number of trials
task.Ntrials = 500;
% number of arrows
task.Nbandits = 4;
% number of observable actions (possible directions)
task.Ndirections = 4;

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'total_iter', 'correct_iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];


%% loop over models
f = {};
f{1} = figure;
f{2} = figure;

nSims = 500;
for model = 0:4
    if model == 4
        %Lazy Arrow
        for it=1:nSims
            data = SimulateNBandits_LA_track_correct(task);
            data = array2table(data);
            data.Properties.VariableNames(1:8+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw = analyzeswitch(data,task);
            rew_data(it,:) = mean(sw.rew);
            cor_data(it, :) = mean(sw.cor);
        end
    end
    
    if model == 3
        %win-stay-lose-shift
        params = [0]; %epsilon (noise)
        % simulate 100 times 
        for it=1:nSims
            data = SimulateNBandits_WSLS_track_correct(task, params);
            data = array2table(data);
            data.Properties.VariableNames(1:8+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw = analyzeswitch(data,task);
            rew_data(it,:) = mean(sw.rew);
            cor_data(it, :) = mean(sw.cor);
        end
    end
    if model ==2
        %RL
        %params = [20,0,.5];
        params = [.5,.2]; %[alpha, stick]
        % simulate 100 times HRL
        for it=1:nSims
            
            data = SimulateNBandits_RL_track_correct(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:8+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw = analyzeswitch(data,task);
            rew_data(it,:) = mean(sw.rew);
            cor_data(it, :) = mean(sw.cor);
        end
        
        
    elseif model==1
        %sticky Bayes
        params = [10,0,2];
        % simulate 100 times HRL - sticky
        for it=1:nSims
            
            data = SimulateNBandits_stick_track_correct(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:8+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw = analyzeswitch(data,task);
            rew_data(it,:) = mean(sw.rew);
            cor_data(it, :) = mean(sw.cor);
        end
    elseif model==0
        %no stick- Bayes
        params = [10,0]; %[beta, epsilon]
        
        % simulate 100 times Inference
        for it=1:nSims
            
            data = SimulateNBandits_track_correct(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:8+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw = analyzeswitch(data,task);
            rew_data(it,:) = mean(sw.rew);
            cor_data(it, :) = mean(sw.cor);
        end
    end
    %plot results
    figure(f{1})
    hold on
    errorbar([-4:9],mean(rew_data),std(rew_data)/sqrt(nSims),'linewidth',2)
    
    figure(f{2})
    hold on
    errorbar([-4:9],mean(cor_data),std(cor_data)/sqrt(nSims),'linewidth',2)
end
figure(f{1})
ylim([0 1])
legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')

figure(f{2})
ylim([0 1])
legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability correct')
xlabel('trials after switch')


%% explore "lazy" heuristic models
figure; hold on
for model = 3:4
  
    if model ==3
        % direction based on # of arrows
        % simulate 100 times 
        for it=1:nSims
            data = SimulateNBandits_LA(task);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    
    if model == 4
        %win-stay-lose-shift side
        params = [0]; %epsilon (noise)
        % simulate 100 times 
        for it=1:nSims
            data = SimulateNBandits_WSLS_side(task, params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    if model == 4
        %win-stay-lose-shift arrow
        params = [0]; %epsilon (noise)
        % simulate 100 times 
        for it=1:nSims
            data = SimulateNBandits_WSLS_arrow(task, params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
end

ylim([.4 1])
legend('Lazy Arrow','WSLS- side', 'WSLS- arrow')
set(gca,'fontsize',14)
%% test effect of parameters
clear all

%% set up task parameters
task = [];

task.prew=[.1 .9];
task.pswitch = .05;
task.Ntrials = 500;
task.Nbandits = 3;

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

%% Focus on HRL
%params = [20,0,.5];


figure;
hold on
% parameter values to explore
ps = [.4:.1:.9]; %vary alpha values (learning rate)
nSims = 100;
% loop over parameter values
for p = ps
    params = [p,.2];
    %params = [p,0];
    % simulate for each parameter value 100 times
    for it=1:nSims
        
        data = SimulateNBandits_RL(task,params);
        data = array2table(data);
        data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
        % analyze and plot behavior
        sw(it,:) = mean(analyzeswitch(data,task));
    end
    
    
    % plot results
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2,'color',p*[0 0 1])
end
ylim([.4 1])
legend('\alpha=.4','.5','.6','.7','.8','.9')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')

%%
nn = 100;
prob = ones(1, nn)/nn;%[0.2, 0.5, 0.3];
opts = 1:nn;
h1 = [];
tic
for i = 1:5000
    h1(end+1) = find(mnrnd(1, prob));
end
t1 = toc;

h2 = [];
tic
for i = 1:5000
    h2(end+1) = randsample(opts, 1, true, prob);
end
t2 = toc;