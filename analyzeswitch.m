function [sw, sw_prob] = analyzeswitch(data,task)
%sw size: number of switches by 14
%gets reward of 4 trials before switch through 9 trials after switch (14
%total)
%also tracks probability of previous cb & new cb

iter = data{:,3};
Nbandits = task.Nbandits;

T = find(iter==1);
T = T(2:end-1); %ignore very first T value since this is very first trial, no switching occured

%r_ind = find(strcmp(data.Properties.VariableNames, 'reward'));

sw = NaN(length(T), 14);
sw_prob = struct;

for b = 1:length(T)
    sw(b,:) = data.reward(T(b)-4:T(b)+9);
    old_cb = data.corr_bandit(T(b)-1);
    new_cb = data.corr_bandit(T(b));
    
    sw_prob.old_bandit(b, :) = data{T(b)-4:T(b)+9, end-Nbandits + old_cb}; %old probability
    sw_prob.new_bandit(b,:) = data{T(b)-4:T(b)+9, end-Nbandits + new_cb}; %new probability
    
end

end