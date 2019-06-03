%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 2
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

%% define time scale
p.dt = 0.1;                 %% Time step (0.1 ms -> 10kHz sampling rate)
ttime = 0:p.dt:425;         %% Duration of the force simulation in time steps

%% define motor neuron pool
p.n = 250;                  %% no of units in pool
p.twtforce_range = 100;
c = 4.191807;               % log3(p.twtforce_range) where 3 is the range of contraction times
p.twtforce = exp( ( log(p.twtforce_range) / (p.n-1)) * [0:p.n-1] );
p.tc = 90 * (1./p.twtforce) .^ (1/c);   %% Calculate time to peak force for all motor units

%% calculate twitch waveforms and 1/2 relaxation time
twt = zeros(length(ttime),p.n);
for i = 1 : p.n
    twt(:,i) = twitch( 0, ttime, p.twtforce(i), p.tc(i) );
    hrt(i) = ttime( min (find ( (twt(:,i)' < (p.twtforce(i)/2)) & (ttime>p.tc(i)) ) ) );
end

% representative motor units to show in twitch plot
reps = [1 60 93 116 138 157 173 186 198 208 217 225 233 239 245 250];

f = figure(2); clf

subplot(221)
hist(p.twtforce, 20);   axis tight
set(gca, 'FontSize', 8);
xlabel('Amplitude (arbitrary units)');
ylabel('Frequency');

a1 = subplot(222);
a2 = copyobj(a1, f);

hist(hrt, 20, 'Parent', a1);
axis(a1, 'tight');
set(a1, 'XAxisLocation', 'top');
set(a1, 'XTick', 100:20:240);
set(get(a1, 'XLabel'), 'String', 'Half Relaxation Time (ms)');
set(a1, 'YTick', []);
set(a1, 'FontSize', 8);

hist(p.tc, 20, 'Parent', a2);
axis(a2, 'tight');
set(a2, 'XTick', 30:10:90);
set(get(a2, 'XLabel'), 'String', 'Contraction Time (ms)');
set(get(a2, 'YLabel'), 'String', 'Frequency');
set(a2, 'FontSize', 8);

subplot(223)
plot(p.twtforce, p.tc, '.');    axis tight
xlabel('Amplitude (arbitrary units)');
ylabel('Contraction Time (ms)');
set(gca, 'FontSize', 8);

a3 = subplot(224);
plot(ttime, twt(:, reps) );
set(a3, 'XLim', [-100 500]);
xlabel('Time (ms)');
ylabel('Twitch Waveform (arbitrary units)');
set(a3, 'FontSize', 8);