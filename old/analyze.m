clear all

%% set up task parameters
task = [];

% probability of reward|incorrect/correct left-right choice
task.prew=[.1 .9];
% probability of correct arrow changing
task.pswitch = .05;
% number of trials
task.Ntrials = 500;
% number of arrows
task.Nbandits = 3;

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};

%% loop over models
figure;
hold on

nSims = 100;
for model = [2, 3, 6]%0:4
    if model == 6
        %RL high alpha
        params = [1,.2]; %[alpha, stick]
        % simulate 100 times HRL
        for it=1:nSims
            
            data = SimulateNBandits_RL(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end

    elseif model == 4
        %Lazy Arrow
        for it=1:nSims
            data = SimulateNBandits_LA(task);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    
    if model == 3
        %win-stay-lose-shift
        params = [0]; %epsilon (noise)
        % simulate 100 times 
        for it=1:nSims
            data = SimulateNBandits_WSLS(task, params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    if model ==2
        %RL
        %params = [20,0,.5];
        params = [.5,.2]; %[alpha, stick]
        % simulate 100 times HRL
        for it=1:nSims
            
            data = SimulateNBandits_RL(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
        
        
    elseif model==1
        %sticky Bayes
        params = [10,0,2];
        % simulate 100 times HRL - sticky
        for it=1:nSims
            
            data = SimulateNBandits_stick(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    elseif model==0
        %no stick- Bayes
        params = [10,0]; %[beta, epsilon]
        
        % simulate 100 times Inference
        for it=1:nSims
            
            data = SimulateNBandits(task,params);
            data = array2table(data);
            data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
            % do data analysis, store it
            sw(it,:) = mean(analyzeswitch(data,task));
        end
    end
    %plot results
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
end
ylim([.4 1])
legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
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