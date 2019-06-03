%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 5
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

% show alternation of N units %

%% define time scale
p.dt = 0.1;                 %% Time step (0.1 ms -> 10kHz sampling rate)
ttime = 0:p.dt:425;         %% Duration of the force simulation in time steps

%% define motor neuron pool
p.n = 250;                  %% no of units in pool
p.twtforce_range = 100;
c = 4.191807;               % log3(p.twtforce_range) where 3 is the range of contraction times
p.twtforce = exp( ( log(p.twtforce_range) / (p.n-1)) * [0:p.n-1] );
p.tc = 90 * (1./p.twtforce) .^ (1/c);   %% Calculate time to peak force for all motor units

N = 3;  % number of units to demo

rand('state',sum(100*clock))                %% randomize random number gen. state

twt = zeros(length(ttime), N);
pu = ceil( rand(N,1) * (p.n - 1) ) + 1;
alt = zeros(length(ttime), 2*N+1);
leg = zeros(2*N+1, 2*N-1);
legchars = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'...
            'R','S','T','U','V','W','X','Y','Z'];

k = 1;

for i = 1 : N
    twt(:,i) = twitch(0, ttime, p.twtforce(pu(i)), p.tc(pu(i)));
end

for n = 1 : N
    c = combnk(1:N,n);
    for i = 1 : size(c,1)
        alt(:,k) = sum(twt(:, c(i,:) ), 2);
        leg(k,1:2:2*n-1) = legchars(c(i,:));
        for m = 2 : 2 : 2 * n - 2
            leg(k,m) = '+';
        end    
        k = k + 1;
    end    
end

figure(5);   clf

plot(ttime,alt,'LineWidth',2);
legend(char(leg),0);
axis tight
xlabel('Time (ms)');
ylabel('Single-Unit and Compound twitchs');