%% Load Data
stim_loc = 'all_stim_22-Feb-2023_9.mat';
mod_loc = 'all_model_data_22-Feb-2023_9.mat';

load(stim_loc)
load(mod_loc)

nSims = stim_params.nSims;
pswitches = stim_params.pswitches;
Ndirections = stim_params.Ndirections;
Nbandits = stim_params.Nbandits;

mod_names = {'Bayes', 'Sticky Bayes', 'RL', 'WSLS', 'Lazy Arrow', 'WSLS-no mem'};
cols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};

%% Plot pstays
figure
bar([mean(raw_model_data(1).pstay_nrew{1, 2})' mean(raw_model_data(2).pstay_nrew{1, 2})']);
ylim([0, 1])
xticks([1:5])
xticklabels({'Bayes', 'StickyBayes', 'RL', 'WSLS', 'LazyArrow'})
legend(strcat('p(switch) = ', string(pswitches)))
title({'p(stay) | no reward', strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']')})

%% Plot curves

plot_vals_pre = NaN(length(betas), length(alphas), length(pswitches)); 
plot_vals_post = NaN(length(betas), length(alphas), length(pswitches));
%%
fig_holder = {};
for p = 1:length(pswitches)
        for f = 1:2
            fig_holder{p, f} = figure;
        end
end
for p = 1:length(pswitches)
    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d);
        
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
 
                for model = 0:4
            
                    % plot results
            figure(fig_holder{p, 1})
            subplot(2, 2, 2*p-1)
            hold on
            errorbar([-4:9],mean(sw),std(sw)/sqrt(nSims),'linewidth',2, 'Color', cols{b})
            mean_sw = mean(sw);
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
    end
figure(fig_holder{1})
ylim([0 1])
legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability of reward')
xlabel('trials after switch')

figure(fig_holder{2})
ylim([0 1])
legend('Bayes','StickyBayes','RL', 'WSLS', 'Lazy Arrow')
set(gca,'fontsize',14)
title({strcat('prew = [', num2str(task.prew(1)), ',', num2str(task.prew(2)), ']'), strcat('pswitch = ', num2str(task.pswitch))})
ylabel('probability correct')
xlabel('trials after switch')

%%
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
                

%%




% loop over parameter values
for p = 1:length(pswitches)
    task.pswitch = pswitches(p);

   for b = 1:length(betas)

        for a = 1:length(alphas)
            
            sw = RL_data_betas(p).reward{b, a};
            sw_prob_old = RL_data_betas(p).prob_old{b, a};
            sw_prob_new = RL_data_betas(p).prob_new{b, a};
            
    
            % plot results
            figure(fig_holder{1, 1})
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

    figure(fig_holder{1, 2})
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
    

    figure(fig_holder{1, 3})
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

%% Contour Plots (Individual)
% get model fits to data
load('model_results_15-Feb-2023_100.mat')

iter = 2; %4 directions, 4 bandits
mods_to_use = 0:3; %exclude Lazy Arrow

contour_figs = {};
for i = 1:3
    contour_figs{i} = figure;
end

for p = 1:length(pswitches)
    for m = 0:length(mods_to_use)- 1
        model = mods_to_use(m+1);
        a_d = linspace(0, 1); %asymptote
        steep = linspace(-10, 10); %steepness
        [x1, y1] = meshgrid(a_d, steep);
        xi = [x1(:) y1(:)];
        
        prob_rew = data_holder(p).rew_coeffs{model+1, iter};
        prob_old_bandit = data_holder(p).old_bandit_coeffs{model+1, iter};
        prob_new_bandit = data_holder(p).new_bandit_coeffs{model+1, iter};

        figure(contour_figs{1})
        subplot(length(mods_to_use), 2, 2*m+p)
        rew_asy = prob_rew(n, 1) + prob_rew(n, 4);
        rew_steep = prob_rew(n, 2);
        [fig_holder,ep]=ksdensity([rew_asy rew_steep],xi);
        X = reshape(ep(:,1),length(a_d),length(steep));
        Y = reshape(ep(:,2),length(a_d),length(steep));
        Z = reshape(fig_holder,length(a_d),length(steep));
        contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
        ylabel('steepness')
        xlabel('asymptote')
        title(mod_names{model+1})

        figure(contour_figs{2})
        subplot(length(mods_to_use), 2, 2*m+p)
        old_asy = prob_old_bandit(n, 1) + prob_old_bandit(n, 4);
        old_steep = prob_old_bandit(n, 2);
        [fig_holder,ep]=ksdensity([old_asy old_steep],xi);
        X = reshape(ep(:,1),length(a_d),length(steep));
        Y = reshape(ep(:,2),length(a_d),length(steep));
        Z = reshape(fig_holder,length(a_d),length(steep));
        contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
        ylabel('steepness')
        xlabel('asymptote')
        title(mod_names{model+1})

        figure(contour_figs{3})
        subplot(length(mods_to_use), 2, 2*m+p)
        new_asy = prob_new_bandit(n, 1) + prob_new_bandit(n, 4);
        new_steep = prob_new_bandit(n, 2);
        [fig_holder,ep]=ksdensity([new_asy new_steep],xi);
        X = reshape(ep(:,1),length(a_d),length(steep));
        Y = reshape(ep(:,2),length(a_d),length(steep));
        Z = reshape(fig_holder,length(a_d),length(steep));
        contour(X,Y,Z,10, 'EdgeColor', cols{model+1},'EdgeAlpha', 0.5)
        ylabel('steepness')
        xlabel('asymptote')
        title(mod_names{model+1})

    end
end

figure(contour_figs{1}); sgtitle('Reward Curve Fit')
figure(contour_figs{2}); sgtitle('Old Bandit Curve Fit')
figure(contour_figs{3}); sgtitle('New Bandit Curve Fit')