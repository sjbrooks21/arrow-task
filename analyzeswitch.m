function [sw_prob, pstay] = analyzeswitch(data,task)
% sw_prob: structure containing arrays with size = number of switches x 14
% gets probability of correct & reward of 4 trials before switch through 9 trials after switch (14
% total)
% also tracks probability of previous correct bandit (prew-switch) & new
% correct bandit (post-switch)
%
% pstay: [nrew, rew]; probability of staying with bandit after no reward/reward

iter = data{:,3};
Nbandits = task.Nbandits;

% pstay
% categorize trials by whether received reward or not
nrew_inds = find(data.reward== 0);
nrew_inds(nrew_inds == task.Ntrials) = []; %remove element if index is last trial
rew_inds = find(data.reward == 1); 
rew_inds(rew_inds == task.Ntrials) = [];

% find frequency stayed with bandit on the following trial
pstay_nrew = sum(data.bandit(nrew_inds + 1) == data.bandit(nrew_inds))/length(nrew_inds);
pstay_rew = sum(data.bandit(rew_inds + 1) == data.bandit(rew_inds))/length(rew_inds);

pstay = [pstay_nrew, pstay_rew];

% find trials where iteration == 1 (switch to new bandit rule)
T = find(iter==1);
T = T(2:end-1); %ignore very first T value since this is very first trial, no switching occured


sw_prob = struct;

for b = 1:length(T)
    old_cb = data.corr_bandit(T(b)-1); %index of old bandit
    new_cb = data.corr_bandit(T(b)); %index of new bandit
    
    sw_prob.correct = data.corr(T(b)-4:T(b)+9);
    sw_prob.reward(b, :) = data.reward(T(b)-4:T(b)+9);
    sw_prob.old_bandit(b, :) = data{T(b)-4:T(b)+9, end-Nbandits + old_cb}; %old probability
    sw_prob.new_bandit(b,:) = data{T(b)-4:T(b)+9, end-Nbandits + new_cb}; %new probability

end

end