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

%stimulus names
stim_names = {};
for i = 1:task.Nbandits
    stim = strcat('stim_', int2str(i));
    stim_names{end+1} = stim;
end

%[t cb iter stim b s cor r]];
col_names = [{'trial', 'corr_bandit', 'iter'}, stim_names]; %, {'bandit', 'side', 'corr', 'reward'}];

%% generate stimuli
stim_params.Nbandits = [3, 4, 5];
stim_params.Ndirections = 4; 
stim_params.pswitches = [0.05, 0.20];
stim_params.weight = 0;
stim_params.nSims = 100;
task.prew = [0.1, 0.9];

prew_str = num2str(task.prew(2));
% uncomment if need to run again
tic
stim = generate_stimuli(task, stim_params);
t = toc; %should take about 8-10 minutes
save(strcat('all_stim_', string(datetime('today')), '_', prew_str(3:end), '.mat'), 'stim_params', 'stim', 'task')

%% run models, save data
%load in stim data, set up params
load('all_stim_2_13_23.mat', 'stim_params', 'stim', 'task')
%load('all_stim_17-Feb-2023_100.mat', 'stim_params', 'stim', 'task')

nSims = stim_params.nSims;
pswitches = stim_params.pswitches;
Ndirections = stim_params.Ndirections;
Nbandits = stim_params.Nbandits;
model_idx = 0:4;
nModels = length(model_idx);

cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
mod_names = {'Bayes', 'Sticky Bayes', 'RL', 'WSLS', 'Lazy Arrow', 'WSLS-no mem'};
mod_params = {[10,0],... %[beta, epsilon]
               [10,0,2],... %[beta, epsilon, stick]
               {[.5, .5], [.2]},... %[alpha, stick]
               [0],... %epsilon (noise)
               [],...%no parameters
               [0]}; %epsilon (noise)

%% Get model data
prew_str = num2str(task.prew(2));
all_model_data = struct;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    all_model_data(p).reward = cell(nModels, length(Ndirections)*length(Nbandits));
    all_model_data(p).prob_old = cell(nModels, length(Ndirections)*length(Nbandits));
    all_model_data(p).prob_new = cell(nModels, length(Ndirections)*length(Nbandits));
    all_model_data(p).pstay_nrew = cell(1, length(Ndirections)*length(Nbandits));
    all_model_data(p).pstay_rew = cell(1, length(Ndirections)*length(Nbandits));
    
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

            all_model_data(p).pstay_nrew{iter} = NaN(nSims, nModels);
            all_model_data(p).pstay_rew{iter}= NaN(nSims, nModels);

            for model = model_idx
                mod_data = get_model_data(task, session_stim, model, mod_params{model+1},col_names);
                all_model_data(p).reward{model+1, iter} = mod_data.rew_data;
                all_model_data(p).prob_old{model+1, iter} = mod_data.prob_old;
                all_model_data(p).prob_new{model+1, iter} = mod_data.prob_new;
                all_model_data(p).pstay_nrew{iter}(:,model+1) = mod_data.pstay(:, 1);
                all_model_data(p).pstay_rew{iter}(:, model+1) = mod_data.pstay(:, 2);
            end
            iter = iter + 1;
        end
    end
end

% mod_pstays = struct;
% 
% for p = 1:length(pswitches)
%     all_mod_pstay_rew = {};
%     for i = 1:5
%         for j = 1:3
%             all_mod_pstay_rew{j}(:, i) = all_model_data(p).pstay{i, j}(:, 2);
%         end
% 
%     end
%     all_mod_pstay_nrew = {};
%     for i = 1:5
%         for j = 1:3
%             all_mod_pstay_nrew{j}(:, i) = all_model_data(p).pstay{i, j}(:, 1);
%         end
%     end
% 
%     mod_pstays(p).pstay_nrew = all_mod_pstay_nrew;
%     mod_pstays(p).pstay_rew = all_mod_pstay_rew;
% 
% end


%save('all_model_data_2_13_23.mat', 'all_model_data')
save(strcat('all_model_data_', string(datetime('today')), '_', prew_str(3:end),'.mat'), 'all_model_data')

%% Plot pstays
figure
hold on

for p = 1:length(pswitches)
    plot(mean(mod_pstays(p).pstay_nrew{1, 2}), '.', 'MarkerSize', 12)
end
    ylim([0, 1])
    xlim([0, 6])
    xticks([1:5])
    xticklabels({'Bayes', 'StickyBayes', 'RL', 'WSLS', 'LazyArrow'})
    legend(string(pswitches))
    title({'p(stay) no reward', strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})


figure
bar([mean(mod_pstays(1).pstay_nrew{1, 2})' mean(mod_pstays(2).pstay_nrew{1, 2})']);
ylim([0, 1])
xticks([1:5])
xticklabels({'Bayes', 'StickyBayes', 'RL', 'WSLS', 'LazyArrow'})
legend(string(pswitches))
title({'p(stay) no reward', strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']')})


%% get model fits to data
%load('all_model_data_2_13_23.mat', 'all_model_data')
%load('all_model_data_17-Feb-2023_1000.mat', 'all_model_data')
data_holder = struct;

logitfittype = fittype('a/(1+exp(-k*(x-b))) + c',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','k', 'b', 'c'});

% options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
% 'Lower',[0 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)
options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
'Lower',[0 -5 -1 0], 'Upper',[2 5 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)

for p = 1:length(pswitches)
    
    task.pswitch = pswitches(p);

    %max(mean over sessions)  
    data_holder(p).max_prew_pre_avg = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).max_prew_post_avg = nan(5, length(Ndirections)*length(Nbandits));

    %max & avg reward of individual sessions
    data_holder(p).prew_pre_ind = cell(5, length(Ndirections)*length(Nbandits));
    data_holder(p).prew_post_ind = cell(5, length(Ndirections)*length(Nbandits));

    %reward fit coeff distribution over session
    data_holder(p).rew_coeffs = cell(5, length(Ndirections)*length(Nbandits));
    data_holder(p).old_bandit_coeffs = cell(5, length(Ndirections)*length(Nbandits));
    data_holder(p).new_bandit_coeffs = cell(5, length(Ndirections)*length(Nbandits));

    %other measures
    data_holder(p).time_forget_old = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).time_learn_new = nan(5, length(Ndirections)*length(Nbandits));
    data_holder(p).pstay = cell(5, length(Ndirections)*length(Nbandits));
        
    iter = 1; %which combination of Nbandits & Ndirections
    
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
  
            for model = 0:4
                sw = all_model_data(p).reward{model+1, iter};
                sw_prob_old = all_model_data(p).prob_old{model+1, iter};
                sw_prob_new = all_model_data(p).prob_new{model+1, iter};
                
                % add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                data_holder(p).max_prew_pre_avg(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                data_holder(p).max_prew_post_avg(model+1, iter) = max(mean_sw(11:end));

                %max and average reward values of each session
                data_holder(p).prew_pre_ind{model+1, iter} = [max(sw(:, 1:4), [], 2) mean(sw(:, 1:4), 2)];
                data_holder(p).prew_post_ind{model+1, iter} = [max(sw(:, 11:end), [], 2) mean(sw(:, 11:end), 2)];

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

                %Curve Fitting 
                
                prob_rew_dist = nan(nSims, 4);
                prob_old_dist = nan(nSims, 4);
                prob_new_dist = nan(nSims, 4);
                

                x = trials;
                x_post = trials(5:end);

                for n = 1:nSims
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

                data_holder(p).rew_coeffs_logit{model+1, iter} = prob_rew_dist;
                data_holder(p).old_bandit_coeffs_logit{model+1, iter} = prob_old_dist;
                data_holder(p).new_bandit_coeffs_logit{model+1, iter} = prob_new_dist;

            end

            iter = iter + 1;
        end
    end

end

save(strcat('model_results_', string(datetime('today')), '_', prew_str(3:end), '.mat'), 'data_holder')
%% Contour Plots (Individual)
iter = 2;
%mods_to_use = [0,1,2]; %Learning
%mods_to_use = [3,4]; %Heuristic
mods_to_use = 0:3;

contour_figs = {};
for i = 1:3
    contour_figs{i} = figure;
end
for p = 1:2
for m = 0:length(mods_to_use)- 1
    model = mods_to_use(m+1);
    a_d = linspace(0, 1); %asymptote
    steep = linspace(-10, 10); %steepness
    [x1, y1] = meshgrid(a_d, steep);
    xi = [x1(:) y1(:)];
    prob_rew_dist = data_holder(p).rew_coeffs_logit{model+1, iter};
    prob_old_dist = data_holder(p).old_bandit_coeffs_logit{model+1, iter};
    prob_new_dist = data_holder(p).new_bandit_coeffs_logit{model+1, iter};

    figure(contour_figs{1})
    subplot(length(mods_to_use), 2, 2*m+p)
    rew_asy = prob_rew_dist(n, 1) + prob_rew_dist(n, 4);
    rew_steep = prob_rew_dist(n, 2);
    [f,ep]=ksdensity([rew_asy rew_steep],xi);
    X = reshape(ep(:,1),length(a_d),length(steep));
    Y = reshape(ep(:,2),length(a_d),length(steep));
    Z = reshape(f,length(a_d),length(steep));
    contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
    ylabel('steepness')
    xlabel('asymptote')
    title(mod_names{model+1})

    figure(contour_figs{2})
    subplot(length(mods_to_use), 2, 2*m+p)
    old_asy = prob_old_dist(n, 1) + prob_old_dist(n, 4);
    old_steep = prob_old_dist(n, 2);
    [f,ep]=ksdensity([old_asy old_steep],xi);
    X = reshape(ep(:,1),length(a_d),length(steep));
    Y = reshape(ep(:,2),length(a_d),length(steep));
    Z = reshape(f,length(a_d),length(steep));
    contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
    ylabel('steepness')
    xlabel('asymptote')
    title(mod_names{model+1})

    figure(contour_figs{3})
    subplot(length(mods_to_use), 2, 2*m+p)
    new_asy = prob_new_dist(n, 1) + prob_new_dist(n, 4);
    new_steep = prob_new_dist(n, 2);
    [f,ep]=ksdensity([new_asy new_steep],xi);
    X = reshape(ep(:,1),length(a_d),length(steep));
    Y = reshape(ep(:,2),length(a_d),length(steep));
    Z = reshape(f,length(a_d),length(steep));
    contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
    ylabel('steepness')
    xlabel('asymptote')
    title(mod_names{model+1})
    
end
end

figure(contour_figs{1}); sgtitle('Reward Curve Fit')
figure(contour_figs{2}); sgtitle('Old Bandit Curve Fit')
figure(contour_figs{3}); sgtitle('New Bandit Curve Fit')


%% RL parameter values to explore
%exploratory so set only do 100 simulations
nSims = 100;

%diff alpha+ and alpha-
neg_vals = [0.30, 0.40, 0.50, 0.60, 0.70];
pos_vals = [0.30, 0.40, 0.50];

task.Nbandits = 4;
task.Ndirections = 4;
iter = 2; %iteration which matches task parameters
model = 2; %looking specifically at RL;

RL_data_t = struct;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);
    for nr = 1:length(neg_vals)
        a_neg = neg_vals(nr);

        for r = 1:length(pos_vals)
            a_pos = pos_vals(r);

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
            session_stim = stim{p}(iter, 1:100);
            rl_params = {[a_neg, a_pos], [.2]}; %learning rates, stickiness

            mod_data = get_model_data(task, session_stim, model, rl_params ,col_names);
            RL_data_t(p).reward{nr, r} = mod_data.rew_data;
            RL_data_t(p).prob_old{nr, r} = mod_data.prob_old;
            RL_data_t(p).prob_new{nr, r} = mod_data.prob_new;

        end
    end
end
%% Plot alpha pairs
rl_fig_holder = {};
for f = 1:2
    rl_fig_holder{1, f} = figure;
end

plot_vals_pre = NaN(length(neg_vals), length(pos_vals), length(pswitches)); 
plot_vals_post = NaN(length(neg_vals), length(pos_vals), length(pswitches));

% loop over parameter values
for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

   for nr = 1:length(neg_vals)

        for r = 1:length(pos_vals)
            
            sw = RL_data_t(p).reward{nr, r};
            sw_prob_old = RL_data_t(p).prob_old{nr, r};
            sw_prob_new = RL_data_t(p).prob_new{nr, r};
            
    
            % plot results
            figure(rl_fig_holder{1, 1})
            subplot(2, 2, 2*p-1)
            hold on
            errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{nr})
            mean_sw = mean(sw);
            plot_vals_pre(nr, r, p) = max(mean_sw(1:4));
            plot_vals_post(nr, r, p) = max(mean_sw(5:end));
            ylabel('probability of reward')
            xlabel('trials after switch')
            ylim([0, 1])

            %probability old/new bandit
            subplot(2, 2, 2*p)
            hold on
            errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{nr})
            errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{nr})
            ylim([0, 1])
            ylabel('probability of bandit')
            xlabel('trials after switch')

        end
    end

    figure(rl_fig_holder{1, 2})
    subplot(2, 2, 2*p-1)
    imagesc(plot_vals_pre(:, :, p));
    set(gca,'fontsize',14)
    title({'pre-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    colormap jet 
    xticks([1 2 3 4 5])
    xticklabels(string(pos_vals))
    xlabel('\alpha_+')
    yticks([1 2 3 4 5])
    yticklabels(string(neg_vals))
    ylabel('\alpha_-')

    subplot(2, 2, 2*p)
    imagesc(plot_vals_post(:, :, p));
    set(gca,'fontsize',14)
    title({'post-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    colormap jet
    xticks([1 2 3 4 5])
    xticklabels(string(pos_vals))
    xlabel('\alpha_+')
    yticks([1 2 3 4 5])
    yticklabels(string(neg_vals))
    ylabel('\alpha_-')
    %set(gca,'fontsize',14)
    %title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
    %legend({'0.25', '0.50', '0.75'})
    %legend({'[0.25, 0.5]', '[0.50, 0.50]', '[0.75, 0.50]'}, 'Location', 'bestoutside')
end
%% RL alpha beta pair vary
betas = [5, 10, 20, 30];
alphas = [0.30, 0.40, 0.50, 0.60, 0.70];

task.Nbandits = 4;
task.Ndirections = 4;
iter = 2; %iteration which matches task parameters
model = 2; %looking specifically at RL;

RL_data_betas = struct;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    for b = 1:length(betas)

        for a = 1:length(alphas)
            alpha = [alphas(a), alphas(a)];

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
            session_stim = stim{p}(iter, 1:100);
            rl_params = {alpha, [0.20], [betas(b)]}; %learning rates, stickiness

            mod_data = get_model_data(task, session_stim, model, rl_params ,col_names);
            RL_data_betas(p).reward{b, a} = mod_data.rew_data;
            RL_data_betas(p).prob_old{b, a} = mod_data.prob_old;
            RL_data_betas(p).prob_new{b, a} = mod_data.prob_new;

        end
    end
end
%% plotting

rl_fig_holder = {};
for f = 1:3
    rl_fig_holder{1, f} = figure;
end

plot_vals_pre = NaN(length(betas), length(alphas), length(pswitches)); 
plot_vals_post = NaN(length(betas), length(alphas), length(pswitches));
beta_colors = linspace(0.1, 1, length(betas));
% loop over parameter values
for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

   for b = 1:length(betas)

        for a = 1:length(alphas)
            
            sw = RL_data_betas(p).reward{b, a};
            sw_prob_old = RL_data_betas(p).prob_old{b, a};
            sw_prob_new = RL_data_betas(p).prob_new{b, a};
            
    
            % plot results
            figure(rl_fig_holder{1, 1})
            subplot(2, 2, 2*p-1)
            hold on
            errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{b})
            mean_sw = mean(sw);
            plot_vals_pre(b, a, p) = max(mean_sw(1:4));
            plot_vals_post(b, a, p) = max(mean_sw(5:end));
            ylabel('probability of reward')
            xlabel('trials after switch')
            ylim([0, 1])

            %probability old/new bandit
            subplot(2, 2, 2*p)
            hold on
            errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{b})
            errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{b})
            ylim([0, 1])
            ylabel('probability of bandit')
            xlabel('trials after switch')

        end
    end

    figure(rl_fig_holder{1, 2})
    subplot(2, 2, 2*p-1)
    imagesc(plot_vals_pre(:, :, p));
    set(gca,'fontsize',14)
    title({'pre-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    colormap parula
    xticks([1 2 3 4 5])
    xticklabels(string(alphas))
    xlabel('\alpha')
    yticks([1 2 3 4])
    yticklabels(string(betas))
    ylabel('\beta')

    subplot(2, 2, 2*p)
    imagesc(plot_vals_post(:, :, p));
    set(gca,'fontsize',14)
    title({'post-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    colormap parula
    xticks([1 2 3 4 5])
    xticklabels(string(alphas))
    xlabel('\alpha')
    yticks([1 2 3 4])
    yticklabels(string(betas))
    ylabel('\beta')
    

    figure(rl_fig_holder{1, 3})
    subplot(2, 2, 2*p-1)
    bar(plot_vals_pre(:, :, p)');
    %plot(alphas, plot_vals_pre(b, :, p), '.', 'MarkerSize', 10, 'Color', cols{b})
    title({'pre-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    xlabel('\alpha')
    xticklabels(string(alphas))
    ylabel('probability of reward')
    ylim([0.5, 1])

    subplot(2, 2, 2*p)
    bar(plot_vals_post(:, :, p)');
    ylim([0.5, 1])
    %plot(alphas, plot_vals_post(b, :, p), '.', 'MarkerSize', 10, 'Color', cols{b})
    title({'post-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    xlabel('\alpha')
    xticklabels(string(alphas))
    ylabel('probability of reward')
    ylim([0.5, 1])
    %legend(string(betas))
    %set(gca,'fontsize',14)
    %title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
    %legend({'0.25', '0.50', '0.75'})
    %legend({'[0.25, 0.5]', '[0.50, 0.50]', '[0.75, 0.50]'}, 'Location', 'bestoutside')
end

%just one example
figure; hold on
b = 3; %beta = 20
colors = linspace(0.1, 1, length(alphas));
for a = 1:length(alphas)
sw = RL_data_betas(p).reward{b, a};
errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', [colors(a), 0, 0])
end
legend(string(alphas))
ylabel('probability of reward')
ylim([0, 1])
xlabel('trials after switch')

figure; hold on
b = 3; %beta = 20
colors = linspace(0.1, 1, length(alphas));
for a = 1:length(alphas)
sw = RL_data_betas(p).reward{b, a};
errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', [colors(a), 0, 0])
end
legend(string(alphas))
ylabel('probability of reward')
ylim([0, 1])
xlabel('trials after switch')

%% Bayes beta sticky pair vary
betas = [5, 10, 20, 30];
stickies = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5];

task.Nbandits = 4;
task.Ndirections = 4;
iter = 2; %iteration which matches task parameters
model = 1; %looking specifically at sticky beta;

Bayes_data_sticky = struct;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    for b = 1:length(betas)

        for s = 1:length(stickies)

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
            session_stim = stim{p}(iter, 1:100);
            bayes_params = [betas(b), 0, stickies(s)]; %beta, epsilon, stickiness

            mod_data = get_model_data(task, session_stim, model, bayes_params ,col_names);
            Bayes_data_sticky(p).reward{b, s} = mod_data.rew_data;
            Bayes_data_sticky(p).prob_old{b, s} = mod_data.prob_old;
            Bayes_data_sticky(p).prob_new{b, s} = mod_data.prob_new;

        end
    end
end
%% plotting

bayes_fig_holder = {};
for f = 1:3
    bayes_fig_holder{1, f} = figure;
end

plot_vals_pre = NaN(length(betas), length(stickies), length(pswitches)); 
plot_vals_post = NaN(length(betas), length(stickies), length(pswitches));

% loop over parameter values
for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

   for b = 1:length(betas)

        for s = 1:length(stickies)
            
            sw = Bayes_data_sticky(p).reward{b, s};
            sw_prob_old = Bayes_data_sticky(p).prob_old{b, s};
            sw_prob_new = Bayes_data_sticky(p).prob_new{b, s};
            
    
            % plot results
            %reward
%             figure(bayes_fig_holder{1, 1})
%             subplot(2, 2, 2*p-1)
%             hold on
%             errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{s})
             mean_sw = mean(sw);
             plot_vals_pre(b, s, p) = max(mean_sw(1:4));
             plot_vals_post(b, s, p) = max(mean_sw(5:end));
%             ylabel('probability of reward')
%             xlabel('trials after switch')
%             ylim([0, 1])

            %probability old/new bandit
%             subplot(2, 2, 2*p)
%             hold on
%             errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2, 'Color', cols{s})
%             errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2, 'Color', cols{s})
%             ylim([0, 1])
%             ylabel('probability of bandit')
%             xlabel('trials after switch')

        end

        figure(bayes_fig_holder{1, 3})
        subplot(2, 2, 2*p-1)
        hold on
        errorbar(stickies, plot_vals_pre(b, :, p), std(plot_vals_pre(b, :, p))/sqrt(nSims), 'linewidth', 2, 'Color', cols{b})
        title({'pre-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
        xlabel('stickiness')
        ylabel('probability of reward')
        ylim([0, 1])

        subplot(2, 2, 2*p)
        hold on
        errorbar(stickies, plot_vals_post(b, :, p), std(plot_vals_post(b, :, p))/sqrt(nSims), 'linewidth', 2, 'Color', cols{b})
        title({'post-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
        xlabel('stickiness')
        ylabel('probability of reward')
        ylim([0, 1])
 
    end

    figure(bayes_fig_holder{1, 2})
    subplot(2, 2, 2*p-1)
    imagesc(plot_vals_pre(:, :, p));
    set(gca,'fontsize',14)
    title({'pre-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    %clim([0 1])
    colormap parula
    xticks(1:length(stickies))
    xticklabels(string(stickies))
    xlabel('stickiness')
    yticks(1:length(betas))
    yticklabels(string(betas))
    ylabel('\beta')

    subplot(2, 2, 2*p)
    imagesc(plot_vals_post(:, :, p));
    set(gca,'fontsize',14)
    title({'post-switch maximum p(reward)', strcat('pswitch = ', num2str(task.pswitch))})
    colorbar
    clim([0.5 1])
    %clim([0 1])
    colormap parula
    xticks(1:length(stickies))
    xticklabels(string(stickies))
    xlabel('stickiness')
    yticks(1:length(betas))
    yticklabels(string(betas))
    ylabel('\beta')

    %set(gca,'fontsize',14)
    %title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
    %legend({'0.25', '0.50', '0.75'})
    %legend({'[0.25, 0.5]', '[0.50, 0.50]', '[0.75, 0.50]'}, 'Location', 'bestoutside')
end


%indiv

b = 2;
        figure
        hold on
        for p = 1:length(pswitches)
            task.pswitch = pswitches(p);
            subplot(1, 2, 1)
            hold on
            plot(stickies, plot_vals_pre(b, :, p), '.-', 'linewidth', 1.2, 'Color', cols{p}, 'MarkerSize', 10)
        title({'pre-switch maximum p(reward)'})
        xlabel('stickiness')
        ylabel('probability of reward')
        ylim([0, 1])

        subplot(1, 2, 2)
        hold on
        plot(stickies, plot_vals_post(b, :, p), '.-', 'linewidth', 1.2, 'Color', cols{p}, 'MarkerSize', 10)
        title({'post-switch maximum p(reward)'})
        xlabel('stickiness')
        ylabel('probability of reward')
        ylim([0, 1])
        end
legend(string(pswitches))
%% WSLS no memory
wsls_model_data = struct;

for p = 1:length(pswitches)
    task.pswitches = pswitches(p);

    wsls_model_data(p).reward = cell(6, length(Ndirections)*length(Nbandits));
    wsls_model_data(p).prob_old = cell(6, length(Ndirections)*length(Nbandits));
    wsls_model_data(p).prob_new = cell(6, length(Ndirections)*length(Nbandits));
    
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

            for model = [3, 5] %just wsls
                mod_data = get_model_data(task, session_stim, model, mod_params{model+1},col_names);
                wsls_model_data(p).reward{model+1, iter} = mod_data.rew_data;
                wsls_model_data(p).prob_old{model+1, iter} = mod_data.prob_old;
                wsls_model_data(p).prob_new{model+1, iter} = mod_data.prob_new;
            end
            iter = iter + 1;
        end
    end
end


wsls_data_holder = struct;

logitfittype = fittype('a/(1+exp(-k*(x-b))) + c',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','k', 'b', 'c'});

% options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
% 'Lower',[0 -10 -1 0], 'Upper',[2 10 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)
options_logit = fitoptions('Method', 'NonlinearLeastSquares', ...
'Lower',[0 -5 -1 0], 'Upper',[2 5 10 1]); %a,k,b,c (a+c = upper asymptote, c = lower asymptote)

for p = 1:length(pswitches)
    
    task.pswitch = pswitches(p);

    %max of the mean  
    wsls_data_holder(p).prew_pre_avg = nan(5, length(Ndirections)*length(Nbandits));
    wsls_data_holder(p).prew_post_avg = nan(5, length(Ndirections)*length(Nbandits));

    %reward (max & avg) distribution over sessions
    wsls_data_holder(p).prew_pre_dist = cell(5, length(Ndirections)*length(Nbandits));
    wsls_data_holder(p).prew_post_dist = cell(5, length(Ndirections)*length(Nbandits));

    %reward fit coeff distribution over session
    wsls_data_holder(p).rew_coeffs_logit = cell(5, length(Ndirections)*length(Nbandits));
    wsls_data_holder(p).old_bandit_coeffs_logit = cell(5, length(Ndirections)*length(Nbandits));
    wsls_data_holder(p).new_bandit_coeffs_logit = cell(5, length(Ndirections)*length(Nbandits));

    %other measures
    wsls_data_holder(p).time_forget_old = nan(5, length(Ndirections)*length(Nbandits));
    wsls_data_holder(p).time_learn_new = nan(5, length(Ndirections)*length(Nbandits));
        
    iter = 1; %which combination of Nbandits & Ndirections
    
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
  
            for model = [3,5]
                sw = wsls_model_data(p).reward{model+1, iter};
                sw_prob_old = wsls_model_data(p).prob_old{model+1, iter};
                sw_prob_new = wsls_model_data(p).prob_new{model+1, iter};
                
                % add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                wsls_data_holder(p).prew_pre_avg(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                wsls_data_holder(p).prew_post_avg(model+1, iter) = max(mean_sw(11:end));

                %max reward distributions
                wsls_data_holder(p).max_prew_pre_dist{model+1, iter} = [max(sw(:, 1:4), [], 2) mean(sw(:, 1:4), 2)];
                wsls_data_holder(p).max_prew_post_dist{model+1, iter} = [max(sw(:, 11:end), [], 2) mean(sw(:, 11:end), 2)];

                %time to forget old bandit
                mean_prob_old = mean(sw_prob_old);
                tfo = find(mean_prob_old(5:end) <= 1/task.Nbandits, 1);
                if any(tfo)
                    wsls_data_holder(p).time_forget_old(model+1, iter) = tfo-1;
                end
                
                %time to learn new bandit
                mean_prob_new = mean(sw_prob_new);
                tln = find(mean_prob_new(5:end) > 0.5, 1);
                if any(tln)
                    wsls_data_holder(p).time_learn_new(model+1, iter) = tln-1;
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

                wsls_data_holder(p).rew_coeffs_logit{model+1, iter} = prob_rew_dist;
                wsls_data_holder(p).old_bandit_coeffs_logit{model+1, iter} = prob_old_dist;
                wsls_data_holder(p).new_bandit_coeffs_logit{model+1, iter} = prob_new_dist;

            end

            iter = iter + 1;
        end
    end

end

%% plot WSLS
fig_holder = {};
n = 1:100;

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

    for f = 1:5 %set up figures
        fig_holder{p, f} = figure;
    end
    
    
    iter = 1; %which combination of Nbandits & Ndirections
    
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
  
            for model = [3, 5]
                model
                sw = wsls_model_data(p).reward{model+1, iter};
                sw_prob_old = wsls_model_data(p).prob_old{model+1, iter};
                sw_prob_new = wsls_model_data(p).prob_new{model+1, iter};

                sw = sw(n, :);
                sw_prob_old = sw_prob_old(n, :);
                sw_prob_new = sw_prob_new(n, :);

                %all together
                figure(fig_holder{p, 1})
                %reward
                subplot(3, 2, 2*iter-1)
                hold on
                errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2)
                ylim([0, 1])
                
                %probability old/new bandit
                subplot(3, 2, 2*iter)
                hold on
                errorbar([-4:9],mean(sw_prob_old),std(sw_prob_old)/sqrt(nSims),'linewidth',2)
                errorbar([-4:9],mean(sw_prob_new),std(sw_prob_new)/sqrt(nSims),'--', 'linewidth',2)
                ylim([0, 1])
                
                % add to data structures
                %max reward probability preswitch
                trials = -4:9;
                mean_sw = mean(sw);
                data_holder(p).prew_pre_avg(model+1, iter) = max(mean_sw(1:4));

                %max reward probability post switch
                data_holder(p).prew_post_avg(model+1, iter) = max(mean_sw(11:end));

                %max reward distributions
                data_holder(p).max_prew_pre_dist{model+1, iter} = [max(sw(:, 1:4), [], 2) mean(sw(:, 1:4), 2)];
                data_holder(p).max_prew_post_dist{model+1, iter} = [max(sw(:, 11:end), [], 2) mean(sw(:, 11:end), 2)];

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

                %Distributions
                prob_rew_dist = wsls_data_holder(p).rew_coeffs_logit{model+1, iter};
                prob_old_dist = wsls_data_holder(p).old_bandit_coeffs_logit{model+1, iter};
                prob_new_dist = wsls_data_holder(p).new_bandit_coeffs_logit{model+1, iter};
                
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
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('Max Probability of Reward')

                %reward fit coeffs
                figure(fig_holder{p, 3})
                subplot(3, 2, 2*iter-1)
                hold on
                pd = fitdist(prob_rew_dist(n, 1) + prob_rew_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('Reward upper asymptote (A+D)')

                subplot(3, 2, 2*iter)
                hold on
                pd = fitdist(prob_rew_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd,linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'LineWidth', 1)
                xlabel('Reward steepness (B)')

                %prob old bandit fit coeffs
                figure(fig_holder{p, 4})
                subplot(3, 3, 3*iter-2)
                hold on
                pd = fitdist(prob_old_dist(n, 1) + prob_old_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('Old Bandit upper asymptote (A+D)')

                subplot(3, 3, 3*iter-1)
                hold on
                pd = fitdist(prob_old_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('Old Bandit lower asymptote (D)')

                subplot(3, 3, 3*iter)
                hold on
                pd = fitdist(prob_old_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'LineWidth', 1)
                xlabel('Old Bandit steepness (B)')

                %prob new bandit fit coeffs
                figure(fig_holder{p, 5})
                subplot(3, 3, 3*iter-2)
                hold on
                pd = fitdist(prob_new_dist(n, 1) + prob_new_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('New Bandit upper asymptote (A+D)')

                subplot(3, 3, 3*iter-1)
                hold on
                pd = fitdist(prob_new_dist(n, 4), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, x_val);
                plot(x_val, rew, 'LineWidth', 1)
                xlabel('New Bandit lower asymptote (D)')

                subplot(3, 3, 3*iter)
                hold on
                pd = fitdist(prob_new_dist(n, 2), 'Kernel');
                %pd = fitdist(max(sw(:, 11:end), [], 2), 'Beta');
                rew = pdf(pd, linspace(-10, 10, 1000));
                plot(linspace(-10, 10, 1000), rew, 'LineWidth', 1)
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
