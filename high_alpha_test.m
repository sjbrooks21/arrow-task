%% Try out RL high alpha
%want to compare WSLS to RL with alpha = 1
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

%% generate stimuli
stim_params.Nbandits = [3, 4, 5];
stim_params.Ndirections = 4; 
stim_params.pswitches = [0.05, 0.20];
stim_params.weight = 0;
stim_params.nSims = 100;
task.prew = [0.1, 0.9];

% uncomment if need to run again
tic
stim = generate_stimuli(task, stim_params);
t = toc; %should take about 8-10 minutes

prew_str = num2str(task.prew(2));
save(strcat('all_stim_high_alpha_', string(datetime('today')), '_', prew_str(3), '.mat'), 'stim_params', 'stim', 'task')

%% set up models to run
nSims = stim_params.nSims;
pswitches = stim_params.pswitches;
Ndirections = stim_params.Ndirections;
Nbandits = stim_params.Nbandits;

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
mod_names = {'Bayes', 'Sticky Bayes', 'RL-normal', 'WSLS', 'Lazy Arrow', 'WSLS-no mem', 'RL-high'};
mod_params = {[10,0],... %[beta, epsilon]
               [10,0,2],... %[beta, epsilon, stick]
               {[0.5, 0.5], [.2]},... %[alpha, stick]
               [0],... %epsilon (noise)
               [],...%no parameters
               [0]...
               {[1, 1], [.2]}}; %epsilon (noise)

%%
RL_high_alpha_model_data = struct;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    RL_high_alpha_model_data(p).reward = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_high_alpha_model_data(p).prob_old = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_high_alpha_model_data(p).prob_new = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    
    iter = 1;

    for d = 1:length(Ndirections)
        task.Ndirections = Ndirections(d);

        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);

            %stimulus names
            stim_names = {};
            for i = 1:task.Nbandits
                stim_name = strcat('stim_', int2str(i));
                stim_names{end+1} = stim_name;
            end

            prob_stim_names = {};
            for i = 1:task.Nbandits
                stim_name = strcat('prob_stim_', int2str(i));
                prob_stim_names{end+1} = stim_name;
            end

            %[t cb iter stim rew b s cor r prob]];
            col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names, {'rew_incorr', 'rew_corr', 'bandit', 'side', 'corr', 'reward'}, prob_stim_names];

            %get stimuli
            session_stim = stim{p}(iter, :);

            for model = [2, 3, 6]
                if model == 6
                    mod_data = get_model_data(task, session_stim, 2, mod_params{model+1},col_names);
                else
                    mod_data = get_model_data(task, session_stim, model, mod_params{model+1},col_names);
                end

                RL_high_alpha_model_data(p).reward{model+1, iter} = mod_data.rew_data;
                RL_high_alpha_model_data(p).prob_old{model+1, iter} = mod_data.prob_old;
                RL_high_alpha_model_data(p).prob_new{model+1, iter} = mod_data.prob_new;
            end
            iter = iter + 1;
        end
    end
end

%save('all_model_data_2_13_23.mat', 'all_model_data')
save(strcat('all_model_data_high_alpha', string(datetime('today')), '_', prew_str(3),'.mat'), 'RL_high_alpha_model_data')

%% get fits- RL high
RL_data_holder_alpha_high = struct;

logitfittype = fittype('a/(1+exp(-k*(x-b))) + c',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','k', 'b', 'c'});

% options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
% 'Lower',[0 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)
options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
'Lower',[0 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)

for p = 1:length(pswitches)
    
    task.pswitch = pswitches(p);

    %max of the mean  
    RL_data_holder_alpha_high(p).prew_pre_avg = nan(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_data_holder_alpha_high(p).prew_post_avg = nan(length(mod_params), length(Ndirections)*length(Nbandits));

    %reward (max & avg) distribution over sessions
    RL_data_holder_alpha_high(p).prew_pre_dist = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_data_holder_alpha_high(p).prew_post_dist = cell(length(mod_params), length(Ndirections)*length(Nbandits));

    %reward fit coeff distribution over session
    RL_data_holder_alpha_high(p).rew_coeffs_logit = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_data_holder_alpha_high(p).old_bandit_coeffs_logit = cell(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_data_holder_alpha_high(p).new_bandit_coeffs_logit = cell(length(mod_params), length(Ndirections)*length(Nbandits));

    %other measures
    RL_data_holder_alpha_high(p).time_forget_old = nan(length(mod_params), length(Ndirections)*length(Nbandits));
    RL_data_holder_alpha_high(p).time_learn_new = nan(length(mod_params), length(Ndirections)*length(Nbandits));
        
    iter = 1; %which combination of Nbandits & Ndirections
    
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
  
            for model = [2, 3, 6]
                sw = RL_high_alpha_model_data(p).reward{model+1, iter};
                sw_prob_old = RL_high_alpha_model_data(p).prob_old{model+1, iter};
                sw_prob_new = RL_high_alpha_model_data(p).prob_new{model+1, iter};
                
                % add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                RL_data_holder_alpha_high(p).prew_pre_avg(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                RL_data_holder_alpha_high(p).prew_post_avg(model+1, iter) = max(mean_sw(11:end));

                %max reward distributions
                RL_data_holder_alpha_high(p).max_prew_pre_dist{model+1, iter} = [max(sw(:, 1:4), [], 2) mean(sw(:, 1:4), 2)];
                RL_data_holder_alpha_high(p).max_prew_post_dist{model+1, iter} = [max(sw(:, 11:end), [], 2) mean(sw(:, 11:end), 2)];

                %time to forget old bandit
                mean_prob_old = mean(sw_prob_old);
                tfo = find(mean_prob_old(5:end) <= 1/task.Nbandits, 1);
                if any(tfo)
                    RL_data_holder_alpha_high(p).time_forget_old(model+1, iter) = tfo-1;
                end
                
                %time to learn new bandit
                mean_prob_new = mean(sw_prob_new);
                tln = find(mean_prob_new(5:end) > 0.5, 1);
                if any(tln)
                    RL_data_holder_alpha_high(p).time_learn_new(model+1, iter) = tln-1;
                end

                %Distributions
                
                prob_rew_dist = nan(nSims, 4);
                prob_old_dist = nan(nSims, 4);
                prob_new_dist = nan(nSims, 4);
                

                x = trials;
                x_post = trials(5:end);

                for n = 1:100%nSims
                    %logistic coeffs
                    

                    y1 = sw(n, 5:end); %reward post-switch
                    fun1_logit = fit(x_post',y1',logitfittype, options_logit);
                    prob_rew_dist(n, :) = coeffvalues(fun1_logit);

                    %y2 = sw_prob_old(n, 5:end); %probability old post-switch
                    y2 = sw_prob_old(n, :); %probability old entire
                    fun2_logit = fit(x',y2',logitfittype, options_logit);
                    prob_old_dist(n, :) = coeffvalues(fun2_logit);

                    %y3 = sw_prob_new(n, 5:end); %probability new post-switch
                    y3 = sw_prob_new(n, :); %probability new entire
                    fun3_logit = fit(x',y3',logitfittype, options_logit);
                    prob_new_dist(n, :) = coeffvalues(fun3_logit);  
                end

                RL_data_holder_alpha_high(p).rew_coeffs_logit{model+1, iter} = prob_rew_dist;
                RL_data_holder_alpha_high(p).old_bandit_coeffs_logit{model+1, iter} = prob_old_dist;
                RL_data_holder_alpha_high(p).new_bandit_coeffs_logit{model+1, iter} = prob_new_dist;

            end

            iter = iter + 1;
        end
    end

end

save(strcat('all_model_data_high_alpha', string(datetime('today')), '_', prew_str(3),'.mat'), 'RL_data_holder_alpha_high', '-append')

%% Plotting

fig_holder = {};
n = 1:nSims;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    for f = 1:5 %set up figures
        fig_holder{p, f} = figure;
    end
    
    
    iter = 1; %which combination of Nbandits & Ndirections
    
    for d = 1 %:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1 %:length(Nbandits)
            task.Nbandits = Nbandits(b);
  
            for model = [2, 3, 6]%0:4
                sw = RL_high_alpha_model_data(p).reward{model+1, iter};
                sw_prob_old = RL_high_alpha_model_data(p).prob_old{model+1, iter};
                sw_prob_new = RL_high_alpha_model_data(p).prob_new{model+1, iter};


                sw = sw(n, :);
                sw_prob_old = sw_prob_old(n, :);
                sw_prob_new = sw_prob_new(n, :);

                %all together
                figure(fig_holder{p, 1})
                %reward
                %subplot(3, 2, 2*iter-1)
                subplot(1, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                
                %probability old/new bandit
                %subplot(3, 2, 2*iter)
                subplot(1, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{model+1})
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{model+1})
                ylim([0, 1])
                
                % add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                RL_data_holder_alpha_high(p).prew_pre_avg(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                RL_data_holder_alpha_high(p).prew_post_avg(model+1, iter) = max(mean_sw(11:end));

                %max reward distributions
                RL_data_holder_alpha_high(p).max_prew_pre_dist{model+1, iter} = [max(sw(:, 1:4), [], 2) mean(sw(:, 1:4), 2)];
                RL_data_holder_alpha_high(p).max_prew_post_dist{model+1, iter} = [max(sw(:, 11:end), [], 2) mean(sw(:, 11:end), 2)];

                %time to forget old bandit
                mean_prob_old = mean(sw_prob_old);
                tfo = find(mean_prob_old(5:end) <= 1/task.Nbandits, 1);
                if any(tfo)
                   RL_data_holder_alpha_high(p).time_forget_old(model+1, iter) = tfo-1;
                end
                
                %time to learn new bandit
                mean_prob_new = mean(sw_prob_new);
                tln = find(mean_prob_new(5:end) > 0.5, 1);
                if any(tln)
                    RL_data_holder_alpha_high(p).time_learn_new(model+1, iter) = tln-1;
                end

                %Distributions
                prob_rew_dist = RL_data_holder_alpha_high(p).rew_coeffs_logit{model+1, iter};
                prob_old_dist = RL_data_holder_alpha_high(p).old_bandit_coeffs_logit{model+1, iter};
                prob_new_dist = RL_data_holder_alpha_high(p).new_bandit_coeffs_logit{model+1, iter};

                
                x = trials;
                x_post = trials(5:end);

                %Distribution plots
                %reward (max post-switch)
                figure(fig_holder{p, 2})
                x_val = linspace(0, 1, 1000);
                subplot(3, 1, iter)
                hold on
                pd = fitdist(max(sw(n, 5:end), [], 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Max Probability of Reward')

                %reward fit coeffs
                figure(fig_holder{p, 3})
                subplot(3, 2, 2*iter-1)
                hold on
                pd = fitdist(prob_rew_dist(n, 1) + prob_rew_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Reward upper asymptote (A+D)')

                subplot(3, 2, 2*iter)
                hold on
                pd = fitdist(prob_rew_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd,linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Reward steepness (B)')

                %prob old bandit fit coeffs
                figure(fig_holder{p, 4})
                subplot(3, 3, 3*iter-2)
                hold on
                pd = fitdist(prob_old_dist(n, 1) + prob_old_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Old Bandit upper asymptote (A+D)')

                subplot(3, 3, 3*iter-1)
                hold on
                pd = fitdist(prob_old_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Old Bandit lower asymptote (D)')

                subplot(3, 3, 3*iter)
                hold on
                pd = fitdist(prob_old_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('Old Bandit steepness (B)')

                %prob new bandit fit coeffs
                figure(fig_holder{p, 5})
                subplot(3, 3, 3*iter-2)
                hold on
                pd = fitdist(prob_new_dist(n, 1) + prob_new_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('New Bandit upper asymptote (A+D)')

                subplot(3, 3, 3*iter-1)
                hold on
                pd = fitdist(prob_new_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('New Bandit lower asymptote (D)')

                subplot(3, 3, 3*iter)
                hold on
                pd = fitdist(prob_new_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'Color',cols{model+1}, 'LineWidth', 1)
                xlabel('New Bandit steepness (B)')

            end

            iter = iter + 1;
        end
    end

    
    for f = 1:5
        figure(fig_holder{p, f})
        sgtitle({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))}) 
    end

end
