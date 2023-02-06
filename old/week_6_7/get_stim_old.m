function data = get_stim(task)
%get task parameters
pswitch = task.pswitch;
Ntrials = task.Ntrials;
Nbandits = task.Nbandits;
Ndirections = task.Ndirections;
prew = task.prew;

data = [];
cb = 1;
iter = 1;


for t = 1:Ntrials
    stim = randi([0,Ndirections-1],[1,Nbandits]);

    rew_incorr = rand<prew(1);
    rew_corr = rand<prew(2);
    
    data = [data;[t cb iter stim rew_incorr rew_corr]];
    
    if iter>10 & rand<pswitch
        iter = 1;
        bs = setdiff(1:Nbandits,cb);
        cb = bs(randi(Nbandits-1));
    else
        iter = iter + 1;
    end
end

%for further analysis
% fun = @(x) numel(unique(x)) == 1;
% d = data(:, 4:6);
% d = num2cell(d, 2);
% all_same = cellfun(fun, d);
% s = zeros(Ntrials, 1);
% s(all_same) = 1;
% figure; 
% h = stem(1:Ntrials, s,'Color','k');
% set(h, 'Marker', 'none')
% title({'Trials with all same direction', 'random'})
% xlabel('trial')
% ylabel('is all same')

end