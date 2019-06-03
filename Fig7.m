%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 7
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

% illustrates the SMUP estimation part of the Incremental Stimulation method

%% select treshold model: diameter-dependent, or direct-depth (random)
%p.ThMod = 'diam';
p.ThMod = 'dirdep';

%% define time scale
p.dt = 0.1;              %% Time step (0.1 ms -> 10kHz sampling rate)
etime = 0:p.dt:10;        %% Duration of the EMG simulation in time steps
ttime = 0:p.dt:425;         %% Duration of the force simulation in time steps
p.del = 5;    % 5 ms delay before stim

p.nsamp = 10;    % number of incremental output changes to record
p.stim_inc = 0.01 / 100;  % incremental stimulation as a fraction of maximum stim

%% define motor neuron pool
p.n = 250;         %% no of units in pool
p.twtforce_range = 100;
c = 4.191807;               % log3(p.twtforce_range) where 3 is the range of contraction times
p.twtforce = exp((log(p.twtforce_range)/(p.n-1))*[0:p.n-1]);
p.tc = 90 * (1./p.twtforce) .^ (1/c);   %% Calculate time to peak force for all motor units
p.rte_cv = 0.0165;    % standard deviation of threshold noise

p.depth_fas = 2 * 1e3;      % fascicle depth (um)
p.diam_min = 5;     p.diam_max = 15;
p.depth_min = p.depth_fas + 200;    p.depth_max = p.depth_fas + 400;

rand('state',sum(100*clock))                %% randomize random number gen. state

% axon diameter spatial distribution in the fascicle
p.depth = normrnd(0,1,1,p.n);   p.depth = (p.depth-min(p.depth))/(max(p.depth)-min(p.depth));
p.depth = p.depth * (p.depth_max - p.depth_min) + p.depth_min;

% axon diameters
p.diam = normrnd(0,1,1,p.n);   p.diam = (p.diam-min(p.diam))/(max(p.diam)-min(p.diam));
p.diam = sort( p.diam * (p.diam_max - p.diam_min) + p.diam_min );

switch p.ThMod
    case 'diam'
        % CALCULATE TRESHOLDS BY DEPTH-DIAMETER MODEL, OR
        disp('Calculating thresholds by depth-diameter model.');
        p.mrte = AxonThreshPairs(p.depth, p.diam);

    case 'dirdep'
        % CALCULATE TRESHOLDS BY DEPTH ONLY (RANDOM)
        disp('Randomizing thresholds.');
        thrmax = AxonThresh(p.depth_max, p.diam_min);
        thrmin = AxonThresh(p.depth_min, p.diam_max);
        p.rrte = p.depth;   p.rrte = (p.rrte-min(p.rrte))/(max(p.rrte)-min(p.rrte));
        p.mrte = p.rrte * (thrmax - thrmin) + thrmin;
end

data = struct('twitch',zeros(length(ttime),p.nsamp),...
              'emg',zeros(length(etime),p.nsamp),...
              'stim',zeros(1,p.nsamp),...
              'nrec',zeros(p.nsamp));

rte_span = max(p.mrte) - min(p.mrte);
minstim = min(p.mrte) * (1 - 3 * p.rte_cv);
maxstim = max(p.mrte) * (1 + 3 * p.rte_cv);
          
%% increment stimulus as each new twitch is identified
j = 0;      rescnt = 0;     stimlev = 0;
%% rescnt means Response Count

disp('Starting Incremental Method.');

while (rescnt < p.nsamp) & (stimlev <= maxstim)

    p.rte = p.mrte + normrnd(0, p.rte_cv*p.mrte);
    stimlev = minstim +  j * p.stim_inc * rte_span;
    %percent_stim = j * p.stim_inc * 100
    [ltw,emg,nrec] = IM_pool(stimlev, p, ttime, etime);       
    
    if nrec
        if ~rescnt
            rescnt = 1
            data.stim(rescnt) = stimlev;
            data.twitch(:,rescnt) = ltw';
            data.emg(:,rescnt) = emg';
            data.nrec(rescnt) = nrec;
        elseif ~ismember( max(ltw), max(data.twitch) )
            rescnt = rescnt + 1
            data.stim(rescnt) = stimlev;
            data.twitch(:,rescnt) = ltw';
            data.emg(:,rescnt) = emg';
            data.nrec(rescnt) = nrec;
        end
    end    

    j = j + 1;

end

if rescnt < p.nsamp
    warning(['All ' num2str(p.nsamp) ' increments were not found - only ' num2str(rescnt) '.']);
    p.nsamp = rescnt;
end

figure(7);  clf
plot(etime, data.emg, 'LineWidth', 2);   axis tight
xlabel('Time (ms)');
ylabel('CMAPs (arbitrary units)');

% to see twitch output of Incremental Stimulation Method...
% plot(ttime, data.twitch, 'LineWidth', 2);