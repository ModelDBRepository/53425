function [tw] = twitch(isi,rtime,force,tc)
%% [tw] = twitch2(isi,rtime,force,tc)
%% generate a twitch, inputs are:
%%   isi = times of spikes
%%   rtime = time sequence for sim
%%   force = peak force for this mu
%%   tc = time constant for this mu
%% output is tw = twitch force over the sim time

tw = zeros(size(rtime));

%% determine time of each spike
stime = isi(isi<rtime(end));

if(length(stime)>0)
    
    [st,rt] = meshgrid(stime,rtime);
    offs = rt-st;
    offs(offs<0)=0;

    %% twitch function for each spike
    twt = ( force*offs.*exp(1-(offs/tc)) )/tc ;
    tw = sum( twt' ,1);   
end
    
tw = tw(:);