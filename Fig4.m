%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 4
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

%% define time scale
p.dt = 0.1;                 %% Time step (0.1 ms -> 10kHz sampling rate)
ttime = 0:p.dt:425;         %% Duration of the force simulation in time steps

%% define motor neuron pool
p.n = 250;                  %% no of units in pool
p.twtforce_range = 100;
p.twtforce = exp( ( log(p.twtforce_range) / (p.n-1)) * [0:p.n-1] );

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

% calculate axon thresholds by model
disp('Calculating thresholds by depth-diameter model.');
p.mrte = AxonThreshPairs(p.depth, p.diam);
        
% calculate axon thresholds by depth (random)
disp('Randomizing thresholds.');
thrmax = AxonThresh(p.depth_max, mean(p.diam));
thrmin = AxonThresh(p.depth_min, mean(p.diam));
p.rrte = p.depth;   p.rrte = (p.rrte-min(p.rrte))/(max(p.rrte)-min(p.rrte));
p.rrte = p.rrte * (thrmax - thrmin) + thrmin;

f = figure(4); clf
set(f, 'Units', 'normalized', 'Position', [0.1 0.55 0.85 0.35]);

subplot(121)
plot(1:p.n, p.mrte, '.', 1:p.n, p.rrte, '.');   axis tight
xlabel('Motor Unit Number');
ylabel('Axon Threshold (mA)');
legend('Diameter (Depolarization Model) Thresholds', 'Depth (Random) Thresholds', 0);

subplot(122)
plot(p.mrte, p.twtforce, '.');    axis tight
xlabel('Axon Threshold (mA)');
ylabel('Twitch Amplitude (arbitrary units)');