function [muap] = muap2_memory(isi,rtime,force)
%% [muap] = muap2(isi,rtime,force)
%% generates a 1st derivative gaussian MUAP
%% saves RAM compared to muap2.m
%% inputs are:
%%   isi = times of spikes
%%   rtime = time sequence for sim
%%   force = peak force for this mu
%%   muap time constant is uniform
%% output is muap = muap over the sim time

muap = zeros(size(rtime));

%% determine time of each spike
stime = isi(isi<rtime(end));

if(length(stime)>0)
        
    for i = 1:length(stime)
        offs = rtime-stime(i);
    
    %% muap function for each spike
    muap = muap - ( sqrt(force)*offs.*exp((-offs.^2)/2) / sqrt(2*pi) );
    
    end

end
    
muap = muap(:);