function [ltw,emg,nrec] = IM_pool(excite,p,ttime,etime)
%%--------------------------------------------------------
%% this function simulates synchronous stimulation of a motor neuron pool 
%% inputs are:
%%       excite = absolute level of excitation
%%       p = params of pool
%%       ttime = time vector for force output
%%       etime = time vector for emg output
%% outputs are:
%%       ltw = muscle force obtained from stimulus
%%       emg = emg output
%%       nrec = number of units recruited for each level of excitation

twlin = zeros(length(ttime),p.n);
emglin = zeros(length(etime),p.n);
nrec = 0;

%% sim for each MN
for i=1:p.n
    
    if excite >= p.rte(i)
        %% do twitch function
        twlin(:,i) = twitch(p.del,ttime,p.twtforce(i),p.tc(i));
        emglin(:,i) = muap2(p.del,etime,p.twtforce(i));
        nrec = nrec + 1;
    end
    
end

ltw = sum(twlin,2)';
emg = sum(emglin,2)';

return