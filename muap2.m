function [muap] = muap2(isi,rtime,force)
%% [muap] = muap2(isi,rtime,force)
%% generate a 1st derivative gaussian MUAP, inputs are:
%%   isi = times of spikes
%%   rtime = time sequence for sim
%%   force = peak force for this mu
%%   muap time constant is uniform
%% output is muap = muap over the sim time

muap = zeros(size(rtime));

%% determine time of each spike
stime = isi(isi<rtime(end));

if(length(stime)>0)
    [st,rt] = meshgrid(stime,rtime);
    offs = rt-st;
    
    %% muap function for each spike
    muapt = -sqrt(force)*offs.*exp((-offs.^2)/2) / sqrt(2*pi);
    
    muap = sum( muapt' ,1);

end
    
muap = muap(:);