function [tw] = twitch_memory(isi,rtime,force,tc)
%% [tw] = twitch_memory(isi,rtime,force,tc)
%% generates muscle twitch using a for loop rather than vectorally
%% saves RAM compared to twitch.m
%% inputs are:
%%   isi = times of spikes
%%   rtime = time sequence for sim
%%   force = peak force for this mu
%%   tc = time constant for this mu
%% output is tw = twitch force over the sim time

tw = zeros(size(rtime));

%% determine time of each spike
stime = isi(isi<rtime(end));

if(length(stime)>0)
    
    for i = 1:length(stime)
        offs = rtime-stime(i);
        offs(offs<0)=0;

    %% twitch function for each spike
        tw = tw + ( force*offs.*exp(1-(offs/tc)) )/tc ;
    
    end
    
end
    
tw = tw(:);