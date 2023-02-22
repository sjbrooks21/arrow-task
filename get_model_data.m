function mod_data = get_model_data(task, session_stim, model, params, col_names)
mod_data = struct;

for it=1:length(session_stim)
    task.stimData = session_stim{it};

    if model == 5 %WSLS no memory
        data = SimulateNBandits_WSLS_dumb(task, params);
    elseif model == 4 %Lazy Arrow
        data = SimulateNBandits_LA(task);
    elseif model == 3 %WSLS
        data = SimulateNBandits_WSLS(task, params);
    elseif model == 2 %RL
        data = SimulateNBandits_RL(task, params);
    elseif model == 1 %Sticky Bayes
        data = SimulateNBandits_stick(task, params);
    elseif model == 0 %Bayes
        data = SimulateNBandits(task, params);
    end

    data = array2table(data);
    data.Properties.VariableNames = col_names;
    % do data analysis, store it
    [rew_data, prob_data] = analyzeswitch(data,task);
    sw(it,:) = mean(rew_data);
    sw_prob_old(it, :) = mean(prob_data.old_bandit);
    sw_prob_new(it, :) = mean(prob_data.new_bandit);
end

mod_data.rew_data = sw;
mod_data.prob_old = sw_prob_old;
mod_data.prob_new = sw_prob_new;

end

