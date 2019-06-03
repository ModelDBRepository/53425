%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 9
%% MODIFIED DATE: 31MAY2005
%%--------------------------------------------------------

% illustrates the SMUP estimation part of the Revised Incremental Stimulation method

% This simulation is time-consuming (2.4h on 2.6GHz Celeron w/ 512MB RAM.
% You can reduce the run time by decreasing the number of SMUP samples
% required (p.nsamp), and/or increasing the stimulus increment (p.stim_inc).

clear

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
p.sweeps = 50;  % number of times to record output at each stimulus level

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

rte_span = max(p.mrte) - min(p.mrte);
minstim = min(p.mrte) * (1 - 4 * p.rte_cv * min(p.mrte));
maxstim = max(p.mrte) * (1 + 3 * p.rte_cv * max(p.mrte));
          
data = struct('twitch',zeros(length(ttime),p.nsamp,p.sweeps),...
              'emg',zeros(length(etime),p.nsamp,p.sweeps),...
              'AFpd',zeros(p.nsamp,1),...
              'AFmx',zeros(p.nsamp,1),...
              'AFmn',zeros(p.nsamp,1),...
              'Fswi',zeros(p.nsamp,1),...
              'APpd',zeros(p.nsamp,1),...
              'APmx',zeros(p.nsamp,1),...
              'APmn',zeros(p.nsamp,1),...
              'Pswi',zeros(p.nsamp,1),...
              'stim',zeros(p.nsamp,1),...
              'nrec',zeros(p.nsamp,p.sweeps));
          
j = 0;      stimlev = 0;		stepcnt = 0;
%% stepcnt means Step Count - the number of observed response steps as
%% defined by the Revised Incremental Stimulation method algorithm

while (stepcnt < p.nsamp) & (stimlev <= maxstim)

    ltw = zeros(p.sweeps,length(ttime));    emg = zeros(p.sweeps,length(etime));
    nrec = zeros(1,p.sweeps);
    
    stimlev = minstim +  j * p.stim_inc * rte_span;
    %percent_stim = j * p.stim_inc * 100

    if ~stepcnt
        for s = 1 : p.sweeps
            p.rte = p.mrte + normrnd(0, p.rte_cv*p.mrte);
            [ltw(s,:),emg(s,:),nrec(s)] = IM_pool(stimlev,p,ttime,etime);
        end    
        if any(nrec)
            stepcnt = 1;
            GetMBBAmps
            if data.AFpd(1) & data.APpd(1) 
                data.stim(stepcnt) = stimlev;
                data.nrec(stepcnt,:) = nrec;
                data.twitch(:,stepcnt,:) = ltw';   data.emg(:,stepcnt,:) = emg';
                stepcnt
            else
                stepcnt = 0;
            end                
        end    
    else
        p.rte = p.mrte + normrnd(0, p.rte_cv*p.mrte);
        [ltw(1,:),emg(1,:),nrec(1)] = IM_pool(stimlev,p,ttime,etime);
        if max(ltw(1,:)) >= data.AFmx(stepcnt);
            for s = 2 : p.sweeps
                p.rte = p.mrte + normrnd(0, p.rte_cv*p.mrte);
                [ltw(s,:),emg(s,:),nrec(s)] = IM_pool(stimlev,p,ttime,etime);
            end		% for s = 2 : p.sweeps
            if all(nrec)
                stepcnt = stepcnt + 1;        
                GetMBBAmps
                if (data.AFmn(stepcnt) >= data.AFmx(stepcnt-1)) & (data.AFpd(stepcnt) > data.AFpd(stepcnt-1))
                    %%% ERROR CATCHING %%%
                        if (data.APmn(stepcnt) < data.APmx(stepcnt-1)) | (data.APpd(stepcnt) <= data.APpd(stepcnt-1))
                            disp('The EMG count-increment doesn''t criteria match the twitch!');
                            disp('Starting this run over.')
                            j = -1;      stimlev = 0;		stepcnt = 0;
                            continue
                        end    
                    %%% ERROR CATCHING %%%
                    data.stim(stepcnt) = stimlev;
                    data.nrec(stepcnt,:) = nrec;
                    data.twitch(:,stepcnt,:) = ltw';   data.emg(:,stepcnt,:) = emg';
                    stepcnt
                else
                    stepcnt = stepcnt - 1;
                    %%% ERROR CATCHING %%%
                        if ~data.AFpd(stepcnt) | ~data.APpd(stepcnt)
                            error('This isn''t the first group.  A predominant response should have been found.');       
                        end     % if ~data.AFpd(stepcnt) | ~data.APpd(stepcnt)
                    %%% ERROR CATCHING %%%
                end         % if (data.AFmn(stepcnt) >= data.AFmx(stepcnt-1)) & (data.AFpd(stepcnt) > data.AFpd(stepcnt-1))
            end             % if all(nrec)    
        end                 % if max(ltw(1,:)) > data.AFmx(stepcnt);
    end                     % if ~stepcnt

    j = j + 1;

end                         % while (stepcnt < p.nsamp) & ((j*p.stim_inc) < 1.2)

%%% ERROR CATCHING %%%
    if ~all(data.AFpd)
        stepcnt
        stimlev
        maxstim
        warning('What happened?  Did the stim level max out?');
    end
%%% ERROR CATCHING %%%

if stepcnt < p.nsamp
    warning(['All ' num2str(p.nsamp) ' increments were not found - only ' num2str(rescnt) '.']);
    p.nsamp = rescnt;
end

twts = zeros(length(ttime),p.nsamp);
pks = zeros(1,p.nsamp);
mups = zeros(length(etime),p.nsamp);
for i = 1 : p.nsamp
    twts(:,i) = data.twitch(:, i, data.Fswi(i));                
    pks(i) = ttime(find(twts(:,i)==data.AFpd(i)));
    mups(:,i) = data.emg(:, i, data.Pswi(i));
end    

figure(9);  clf
for s = 1 : p.sweeps
    for r = 1 : 2 : p.nsamp
        plot(etime,data.emg(:,r,s), 'Color', 0.75*ones(1,3), 'LineWidth', 0.5); hold on
    end    
    for r = 2 : 2 : p.nsamp
        plot(etime,data.emg(:,r,s), 'Color', 0.5*ones(1,3), 'LineWidth', 0.5); hold on
    end    
end
plot(etime,mups,'k','LineWidth',1.5); hold off
xlabel('Time (ms)');
ylabel('CMAPs (arbitrary units)');

% to see twitch output of Revised Incremental Stimulation Method...
% replace etime with ttime, replace .emg with .twitch, repace mups with twts