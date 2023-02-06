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
% number of observable actions (possible directions)
task.Ndirections = 2;

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

%% loop over models
Nbandits = [3, 4, 5];
Ndirections = 4; %[2, 3, 4];

pswitches = [0.05, 0.15, 0.20, 0.3]; %[0.05, 0.10, 0.15, 0.20];

nSims = 100;

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
mod_names = {'Bayes', 'Sticky Bayes', 'RL', 'WSLS', 'Lazy Arrow', 'All'};
fig_holder = {};
data_holder = {};

for p = 1:length(pswitches)
    
    data_holder(p) = struct;
    
    for f = 1:6 %set up figuresS
        fig_holder{p, f} = figure;
    end
    
    %r_plot = figure;
    %p_plot = figure;
    
    task.pswitch = pswitches(p);
    
    iter = 1; %which combination of Nbandits & Ndirections
    for b = 1:length(Nbandits)
        task.Nbandits = Nbandits(b);
        
        for d = 1:length(Ndirections)
            
            task.Ndirections = Ndirections(d);
            
            %stimulus names
            stim_names = {};
            for i = 1:task.Nbandits
                stim = strcat('stim_', int2str(i));
                stim_names{end+1} = stim;
            end
            
            prob_stim_names = {};
            for i = 1:task.Nbandits
                stim = strcat('prob_stim_', int2str(i));
                prob_stim_names{end+1} = stim;
            end
            
            %[t cb iter stim rew b s cor r prob]];
            col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names, {'rew_incorr', 'rew_corr', 'bandit', 'side', 'corr', 'reward'}, prob_stim_names];
            
            %get stimuli
            session_stim = {};
            for i = 1:nSims
                session_stim{i} = get_stim(task);
            end
            
            for model = 0:4
                if model == 4
                    %Lazy Arrow
                    for it=1:nSims
                        task.stimData = session_stim{it};
                        data = SimulateNBandits_LA(task);
                        data = array2table(data);
                        data.Properties.VariableNames = col_names;
                        % do data analysis, store it
                        [rew_data, prob_data] = analyzeswitch(data,task);
                        sw(it,:) = mean(rew_data);
                        sw_prob_old(it, :) = mean(prob_data.old_bandit);
                        sw_prob_new(it, :) = mean(prob_data.new_bandit);
                    end
                end
                
                if model == 3
                    %win-stay-lose-shift
                    params = [0]; %epsilon (noise)
                    % simulate 100 times
                    for it=1:nSims
                        task.stimData = session_stim{it};
                        data = SimulateNBandits_WSLS(task, params);
                        data = array2table(data);
                        data.Properties.VariableNames = col_names;
                        % do data analysis, store it
                        [rew_data, prob_data] = analyzeswitch(data,task);
                        sw(it,:) = mean(rew_data);
                        sw_prob_old(it, :) = mean(prob_data.old_bandit);
                        sw_prob_new(it, :) = mean(prob_data.new_bandit);
                    end
                end
                if model ==2
                    %RL
                    params = [.5,.2]; %[alpha, stick]
                    % simulate 100 times HRL
                    for it=1:nSims
                        task.stimData = session_stim{it};
                        data = SimulateNBandits_RL(task,params);
                        data = array2table(data);
                        data.Properties.VariableNames = col_names;
                        % do data analysis, store it
                        [rew_data, prob_data] = analyzeswitch(data,task);
                        sw(it,:) = mean(rew_data);
                        sw_prob_old(it, :) = mean(prob_data.old_bandit);
                        sw_prob_new(it, :) = mean(prob_data.new_bandit);
                    end
                    
                    
                elseif model==1
                    %sticky Bayes
                    params = [10,0,2];
                    for it=1:nSims
                        task.stimData = session_stim{it};
                        data = SimulateNBandits_stick(task,params);
                        data = array2table(data);
                        data.Properties.VariableNames = col_names;
                        % do data analysis, store it
                        [rew_data, prob_data] = analyzeswitch(data,task);
                        sw(it,:) = mean(rew_data);
                        sw_prob_old(it, :) = mean(prob_data.old_bandit);
                        sw_prob_new(it, :) = mean(prob_data.new_bandit);
                    end
                elseif model==0
                    %no stick- Bayes
                    params = [10,0]; %[beta, epsilon]
                    
                    % simulate 100 times Inference
                    for it=1:nSims
                        task.stimData = session_stim{it};
                        data = SimulateNBandits(task,params);
                        data = array2table(data);
                        data.Properties.VariableNames = col_names;
                        % do data analysis, store it
                        [rew_data, prob_data] = analyzeswitch(data,task);
                        sw(it,:) = mean(rew_data);
                        sw_prob_old(it, :) = mean(prob_data.old_bandit);
                        sw_prob_new(it, :) = mean(prob_data.new_bandit);
                    end
                end
                %plot results
                %errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
                %                 figure(r_plot)
                %                 subplot(3, 1, iter)
                %                 hold on
                %                 errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
                %
                %                 figure(p_plot)
                %                 subplot(3, 1, iter)
                %                 hold on
                %                 errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                %                 errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{model+1})
                
                %specific model
                figure(fig_holder{p, model+1})
                %reward
                subplot(3, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                %probability
                subplot(3, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                
                %all together
                figure(fig_holder{p, 6})
                %reward
                subplot(3, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                
                %probability
                subplot(3, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
            end
            %ylim([.20 1])
            %legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
            %set(gca,'fontsize',14)
            %title({strcat('Nbandits = ', num2str(task.Nbandits)), strcat('Ndirections = ', num2str(task.Ndirections))})
            
            
            %             figure(r_plot)
            %             ylabel('probability of reward')
            %             xlabel('trials after switch')
            %
            %             figure(p_plot)
            %             ylabel('probability of bandit')
            %             xlabel('trials after switch')
            
            iter = iter + 1;
        end
    end
    
    for f = 1:6
        figure(fig_holder{p, f})
        sgtitle({mod_names{f}, strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
    end
end
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
    
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
end

ylim([.4 1])
legend('Lazy Arrow','WSLS')
set(gca,'fontsize',14)
%% test effect of model parameters
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

%% test effect of model parameters
clear all

%% set up task parameters
task = [];

task.prew=[.1 .9];
task.pswitch = .05;
task.Ntrials = 500;
task.Nbandits = 4;

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

%% Focus on HRL- alpha
% %params = [20,0,.5];
%
%
% figure;
% hold on
% % parameter values to explore
% ps = [.4:.1:.9]; %vary alpha values (learning rate)
% nSims = 100;
% % loop over parameter values
% for p = ps
%     params = [p,.2];
%     %params = [p,0];
%     % simulate for each parameter value 100 times
%     for it=1:nSims
%
%         data = SimulateNBandits_RL(task,params);
%         data = array2table(data);
%         data.Properties.VariableNames(1:7+task.Nbandits) = [col_names, {'bandit', 'side', 'corr', 'reward'}];
%         % analyze and plot behavior
%         sw(it,:) = mean(analyzeswitch(data,task));
%     end
%
%
%     % plot results
%     errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2,'color',p*[0 0 1])
% end
% ylim([.4 1])
% legend('\alpha=.4','.5','.6','.7','.8','.9')
% set(gca,'fontsize',14)
% title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
% ylabel('probability of reward')
% xlabel('trials after switch')

%% Focus on HRL- beta
%params = [20,0,.5];


figure;
hold on
% parameter values to explore
ps = [1, 2, 5, 10, 20, 30]; %vary beta values (stochasicity)
colors = linspace(0.1, 1, length(ps));
ii = 1;
nSims = 100;
% loop over parameter values
for p = ps
    params = [0.5, .2, p];
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
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'color',[colors(ii) 0 0])
    ii = ii + 1;
end
ylim([.4 1])
legend('\beta= 1', '2', '5', '10','20','30')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')
%% Focus on HRL- epsilon
%params = [20,0,.5];


figure;
hold on
% parameter values to explore
ps = [0.05:0.1:0.55]; %vary beta values (stochasicity)
nSims = 100;
% loop over parameter values
for p = ps
    params = [0.5, .2, p];
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
legend('\epsilon= 0.05','0.15','0.25','0.35','0.45', '0.55')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')

%% Focus on HRL- beta & vary direction, bandits
Nbandits = [3, 4, 5];
Ndirections = [2, 3, 4];

pswitches = [0.05, 0.20]; %[0.05, 0.10, 0.15, 0.20];
ps = [1, 2, 5, 10, 20, 30]; %vary beta values (stochasicity)
colors = linspace(0.1, 1, length(ps));

for pswitch = 1:length(pswitches)
    
    task.pswitch = pswitches(pswitch);
    
    figure
    
    iter = 1;
    for b = 1:length(Nbandits)
        task.Nbandits = Nbandits(b);
        
        for d = 1:length(Ndirections)
            
            task.Ndirections = Ndirections(d);
            
            %stimulus names
            stim_names = {};
            for i = 1:task.Nbandits
                stim = strcat('stim_', int2str(i));
                stim_names{end+1} = stim;
            end
            
            %[t cb iter stim b s cor r]];
            col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];
            
            nSims = 100;
            
            subplot(3, 3, iter)
            hold on
            
            for model = 2%0:4
                if model == 4
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
                    ii = 1;
                    nSims = 100;
                    % loop over parameter values
                    for p = ps
                        params = [0.5, .2, p];
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
                        errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'color',[colors(ii) 0 0])
                        ii = ii + 1;
                    end
                    ylim([.4 1])
                    %legend('\beta= 1', '2', '5', '10','20','30')
                    set(gca,'fontsize',14)
                    %title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
                    ylabel('probability of reward')
                    xlabel('trials after switch')
                end
                
                %plot results
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
            end
        end
    end
    
    sgtitle({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
end

%% parameter values to explore
ps = [1, 2, 5, 10, 20, 30]; %vary beta values (stochasicity)
colors = linspace(0.1, 1, length(ps));
ii = 1;
nSims = 100;
% loop over parameter values
for p = ps
    params = [0.5, .2, p];
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
    errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'color',[colors(ii) 0 0])
    ii = ii + 1;
end
ylim([.4 1])
%legend('\beta= 1', '2', '5', '10','20','30')
set(gca,'fontsize',14)
%title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')