function sw = analyzeswitch(data,task)
%sw size: number of switches by 14
%gets reward of 4 trials before switch through 9 trials after switch (14
%total)

iter = data{:,3};

T = find(iter==1);
T = T(2:end-1); %ignore very first T value since this is very first trial, no switching occured

r_ind = find(strcmp(data.Properties.VariableNames, 'reward'));

for b = 1:length(T)
    sw(b,:) = data{T(b)-4:T(b)+9, r_ind};
end

end