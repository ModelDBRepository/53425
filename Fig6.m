%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 6
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

skip = 20;    % number of units to skip in between probability plots

%% select treshold model: diameter-dependent, or direct-depth (random)
%p.ThMod = 'diam';
p.ThMod = 'dirdep';

%% define time scale
p.n = 250;         %% no of units in pool
p.rte_cv = 0.0165;    % standard deviation of threshold noise

%% define motor neuron pool
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

ucdf = unique([[1 : skip: p.n] p.n]);      % select units to plot
ncdfs = length(ucdf);       % number of units in the pool to plot

switch p.ThMod
    case 'diam'
        % CALCULATE TRESHOLDS BY DEPTH-DIAMETER MODEL, OR
        disp('Calculating thresholds by depth-diameter model.');
        p.mrte = sort( AxonThreshPairs(p.depth, p.diam) );
        dmrte = p.mrte(ucdf);
    
    case 'dirdep'
        % CALCULATE TRESHOLDS BY DEPTH ONLY (RANDOM)
        disp('Randomizing thresholds.');
        thrmax = AxonThresh(p.depth_max, p.diam_min);
        thrmin = AxonThresh(p.depth_min, p.diam_max);
        p.rrte = p.depth;   p.rrte = (p.rrte-min(p.rrte))/(max(p.rrte)-min(p.rrte));
        p.rrte = sort(p.rrte) * (thrmax - thrmin) + thrmin;
        dmrte = p.rrte(ucdf);
end

minstim = min(dmrte) * (1 - 3 * p.rte_cv);
maxd = max(dmrte) * (1 + 3 * p.rte_cv);
rr = linspace(minstim, maxd, 5000);
alt_cdf = zeros(ncdfs, length(rr));

% calculate cummulative probability distributino for the selected units
for i = 1 : ncdfs
    alt_cdf(i,:) = cdf('norm', rr, dmrte(i), p.rte_cv*dmrte(i) );
end

figure(6); clf

plot(rr, alt_cdf', 'LineWidth', 2);     axis tight;
xlabel('Axon Threshold (mA)');
ylabel('Probability of Activation');
legend(num2str(ucdf'), 0);