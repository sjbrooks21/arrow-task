function all_stims = generate_stimuli(task, stim_params)

Nbandits = stim_params.Nbandits;
Ndirections = stim_params.Ndirections;
pswitches = stim_params.pswitches;
nSims = stim_params.nSims;
all_stims = {};

for p = 1:length(pswitches)
    task.pswitch = pswitches(p);
    session_stim = {};
    iter = 1;

    for d = 1:length(Ndirections) 
        task.Ndirections = Ndirections(d); 
        for b = 1:length(Nbandits)
            task.Nbandits = Nbandits(b);
            for i = 1:nSims
                session_stim{iter, i} = get_stim3(task);
            end 
            iter = iter + 1;
        end
    end
    all_stims{p} = session_stim;
end

