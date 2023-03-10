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
Nbandits = 4;%[3, 4, 5];
Ndirections = 4; %[2, 3, 4];

pswitches = [0.05]; %, 0.15, 0.20, 0.3]; %[0.05, 0.10, 0.15, 0.20];

%ratio of all same:all unique: some different
r_sud = {[1, 1, 1], [1, 6, 57], [0, 1, 0]}; 

%decide whether to use weights or not
weight = 1;

nSims = 100;

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
mod_names = {'Bayes', 'Sticky Bayes', 'RL', 'WSLS', 'Lazy Arrow', 'All'};
fig_holder = {};

data_holder = struct;

logitfittype = fittype('a/(1+exp(-k*(x-b))) + c',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','k', 'b', 'c'});

options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
'Lower',[-2 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)

%for p = 1:length(pswitches)
    
    data_holder(p).prew_pre = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).prew_post = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).time_forget_old = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).time_learn_new = nan(5, length(Ndirections)*length(Nbandits));

    data_holder(p).rew_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));
    data_holder(p).old_bandit_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));
    data_holder(p).new_bandit_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));

    for f = 1:5+length(r_sud) %set up figures
        fig_holder{p, f} = figure;
    end
    
    task.pswitch = pswitches(p);
    
    for r = 1:length(r_sud)
    iter = 1; %which combination of Nbandits & Ndirections
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);

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
                session_stim{i} = get_stim2(task, weight, r_sud{r}, 0);
                %session_stim{i} = get_stim2(task);
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
                    params = {[0.5,0.5], [0.2]}; %[alpha, stick]
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
                    
                    % simulate 100 times 
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
                
                %specific model
                figure(fig_holder{p, model+1})
                %reward
                subplot(3, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{r})
                ylim([0, 1])
                %probability
                subplot(3, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{r})
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{r})
                ylim([0, 1])
                sgtitle(mod_names{model + 1})

                %all together
                %figure(fig_holder{p, 6})
                figure(fig_holder{p, 5 + r})
                %reward
                subplot(3, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                
                %probability old/new bandit
                subplot(3, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])

                sgtitle(strcat('sdr = ', num2str(r_sud{r})))
                
                %add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                data_holder(p).prew(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                data_holder(p).prew(model+1, iter) = max(mean_sw(1:4));

                %time to forget old bandit
                mean_prob_old = mean(sw_prob_old);
                tfo = find(mean_prob_old(5:end) <= 1/task.Nbandits, 1);
                if any(tfo)
                    data_holder(p).time_forget_old(model+1, iter) = tfo-1;
                end
                
                %time to learn new bandit
                mean_prob_new = mean(sw_prob_new);
                tln = find(mean_prob_new(5:end) > 0.5, 1);
                if any(tln)
                    data_holder(p).time_learn_new(model+1, iter) = tln-1;
                end

%                  %exponential/logistic coeffs
%                  x = trials(5:end);
% 
%                  y1 = mean_sw(5:end); %reward post-switch
%                  fun1_logit = fit(x',y1',logitfittype, options_logit);
%                  data_holder(p).rew_coeffs_logit(model+1, :, iter) = coeffvalues(fun1_logit);
%                  %figure; plot(fun1_logit, x, y1)
%  
%                  y2 = mean_prob_old(5:end); %probability old post-switch 
%                  fun2_logit = fit(x',y2',logitfittype, options_logit);
%                  data_holder(p).old_bandit_coeffs_logit(model+1, :, iter) = coeffvalues(fun2_logit);
%                  %figure; plot(fun2_logit, x, y2)
% 
%                  y3 = mean_prob_new(5:end); %probability new post-switch 
%                  fun3_logit = fit(x',y3',logitfittype, options_logit);
%                  data_holder(p).new_bandit_coeffs_logit(model+1, :, iter) = coeffvalues(fun3_logit);
    
            end
            
            
            iter = iter + 1;
        end
    end

    
%     figure(fig_holder{p, 5 + r})
%     if weight
%         sgtitle(strcat('sdr = ', num2str(r_sud{r})));
%     else
%         sgtitle("No weighting");
%     end
    %sgtitle({'Exponential', strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('sdr = ', num2str(sdr(p)))});%, strcat('pswitch = ', num2str(task.pswitch))})

    end
%end

%% just test sticky bayes over p switches
Nbandits = [3, 4];%[3, 4, 5];
Ndirections = 4; %[2, 3, 4];

pswitches = [0.05, 0.10, 0.15, 0.20, 0.25]; %, 0.15, 0.20, 0.3]; %[0.05, 0.10, 0.15, 0.20];

colors = linspace(0.1, 1, length(pswitches));

%ratio of all same:all unique: some different
r_sud = {[1, 6, 57]}; 

%decide whether to use weights or not
weight = 0;

nSims = 100;

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
mod_names = {'Bayes', 'Sticky Bayes', 'RL', 'WSLS', 'Lazy Arrow', 'All'};
fig_holder = {};

data_holder = struct;

logitfittype = fittype('a/(1+exp(-k*(x-b))) + c',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','k', 'b', 'c'});

options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
'Lower',[-2 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)

figure

for p = 1:length(pswitches)
    
    data_holder(p).prew_pre = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).prew_post = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).time_forget_old = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).time_learn_new = nan(5, length(Ndirections)*length(Nbandits));

    data_holder(p).rew_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));
    data_holder(p).old_bandit_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));
    data_holder(p).new_bandit_coeffs_logit = nan(5, 4, length(Ndirections)*length(Nbandits));

    
    task.pswitch = pswitches(p);
    

    iter = 1; %which combination of Nbandits & Ndirections
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);

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
                session_stim{i} = get_stim2(task, weight, r_sud{r}, 0);
                %session_stim{i} = get_stim2(task);
            end
            
            for model = 1
                
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
                    params = {[.5,0.5], [0.2]}; %[alpha, stick]
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
                    
                    % simulate 100 times 
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
                
                %reward
                subplot(2, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', [0 colors(p) 0])
                ylim([0, 1])
                %probability
                subplot(2, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', [0 colors(p) 0])
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', [0 colors(p) 0])
                ylim([0, 1])
                sgtitle(mod_names{model + 1})

                %add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                data_holder(p).prew(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                data_holder(p).prew(model+1, iter) = max(mean_sw(1:4));

                %time to forget old bandit
                mean_prob_old = mean(sw_prob_old);
                tfo = find(mean_prob_old(5:end) <= 1/task.Nbandits, 1);
                if any(tfo)
                    data_holder(p).time_forget_old(model+1, iter) = tfo-1;
                end
                
                %time to learn new bandit
                mean_prob_new = mean(sw_prob_new);
                tln = find(mean_prob_new(5:end) > 0.5, 1);
                if any(tln)
                    data_holder(p).time_learn_new(model+1, iter) = tln-1;
                end
    
            end
            
            
            iter = iter + 1;
        end
    end

    
%     figure(fig_holder{p, 5 + r})
%     if weight
%         sgtitle(strcat('sdr = ', num2str(r_sud{r})));
%     else
%         sgtitle("No weighting");
%     end
    %sgtitle({'Exponential', strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('sdr = ', num2str(sdr(p)))});%, strcat('pswitch = ', num2str(task.pswitch))})

    
end