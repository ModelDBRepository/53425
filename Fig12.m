%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 12
%% MODIFIED DATE: 31MAY2005
%%--------------------------------------------------------

% illustrates the SMUP estimation part of the Spike-Triggered Averaging method

% This simulation is time-consuming (XXX on 2.6GHz Celeron w/ 512MB RAM.
% You can reduce the run time by decreasing the number of SMUP samples
% required (p.nsta), and/or decreasing the duration of each contraction
% (rtime) which will result in fewer triggers per average.

%% If OUT OF MEMORY during function call STA_pool2.m, replace
%% twitch.m with twitch_memory.m (within STA_pool2.m), and replace
%% muap2.m with muap2_memory.m (within STA_pool2.m).
%% This alternative will increase the simulation run time to about XXX.

clear

%% define time scale
p.dt = 0.5;              %% Time step (0.5 ms -> 2kHz sampling rate)
rtime = 0:p.dt:25000;       %% Duration of the whole simulation in time steps
ttime = -5:p.dt:425;       %% Duration of the force simulation in time steps
etime = 0:p.dt:10;        %% Duration of the EMG simulation in time steps

p.nsta = 20;        % number of motor units to sample
p.maxperMVC = 10 / 100;     % maximum fraction of MVC
p.frmax = 10;       % maximum firing rate of units to sample
p.del = 5;          % 5 ms delay before 100% stim

%% define motor neuron pool
p.n = 250;         %% no of units in pool
p.twtforce_range = 100;  %% force range 
p.twtforce = exp((log(p.twtforce_range)/(p.n-1))*[0:p.n-1]);  %% Calculate peak force for all motor units
c = 4.191807;      % log3(p.force_range) where 3 is the range of contraction times
p.tc = 90*(1./p.twtforce).^(1/c);   %% Calculate time to peak force for all motor units

p.rte_range = 30;  %% recruitment range 
p.rte = exp((log(p.rte_range)/(p.n-1))*[0:p.n-1]); %% Calculate recruitment for all motor units

%% define firing rates
p.mfr = 8;    %% min firing rate
p.gain = 1.5;  %% gain on firing rate
p.pfr = 45-10*p.rte/p.rte(end);   %% peak firing rate (35 down to 45 Hz)
p.cv = 0.2;   %% coeff of variation (use p.cv = 0 for synchronous firing)

p.maxe = p.rte(end)+(p.pfr(end)-p.mfr)/p.gain;  %% max excitation

MVCquicktest                                        %% generate a lookup table of excitation versus force
p.maxlev = interp1(mF/max(mF),lev,p.maxperMVC);      %% interpolate to find excitation for specified %MVC
disp(['Excitation required to achieve ' num2str(p.maxperMVC*100) '% MVC is ' num2str(round(p.maxlev*100)) '% of maximum.']);

rand('state',sum(100*clock))                        %% randomize random number gen. state

p.lev = rand(p.nsta, 1) * p.maxlev;                 %% levels of excitation to use (fraction of p.maxe)
p.lev(p.lev < (p.rte(1)/p.maxe)) = p.rte(1)/p.maxe;

ltw = zeros(p.nsta,length(rtime));    emg = ltw;
p.sMUs = zeros(p.nsta,1);     mfrate = p.sMUs;
numspk = round(rtime(end)/1000*p.pfr(1)*1.05);      %% maximum number of spikes in spike train
ntrigs = round(rtime(end)/1000*p.frmax*1.05);       %% maximum number of triggers in spike train
trigs = zeros(p.nsta,ntrigs);

for j = 1 : p.nsta

    while ~p.sMUs(j)    
        contraction_number = j
        percent_excitation = p.lev(j)*100
        [ltw(j,:),emg(j,:),p.sMUs(j),mfrate(j),trigs(j,:)] = STA_pool2(p.lev(j), p, rtime, numspk, ntrigs);      
        if ~p.sMUs(j) 
        	p.lev(j) = 0;
            while p.lev(j) < (p.rte(1)/p.maxe)
			    p.lev(j) = rand * p.maxlev;    
            end
        end
    end     
        
end
        
trigs(trigs>rtime(end)) = 0;

cutoff = 0.5;
filtorder = 5;
smplrate = 1000/p.dt;
[b,a]=butter(filtorder,2*cutoff/smplrate,'high');

ss = 300/p.dt:length(rtime);

Fw = length(ttime)-1;                % maximum twitch duration
Pw = length(etime)-1;                % muap duration is about 8 ms

Fwin = zeros(Fw+1,size(trigs,2));   Fsmu = zeros(Fw+1,p.nsta);
Pwin = zeros(Pw+1,size(trigs,2));   Psmu = zeros(Pw+1,p.nsta);

Fcount = zeros(p.nsta,1);           Pcount = Fcount;

for j = 1 : p.nsta
    
    disp(['Calculating averages for STA number ' num2str(j) ' of ' num2str(p.nsta) '.']);
    Fdone = 0;  Pdone = 0;
    
    fltw = filtfilt(b,a,ltw(j,ss:end));

    Ftrgs = trigs(j, trigs(j,:) > (ss(1)*p.dt-ttime(1))) + ttime(1);
    Ptrgs = trigs(j, trigs(j,:) > (Pw/2*p.dt) );
    
    for s = 1 : max([length(Ftrgs) length(Ptrgs)])
     
        if ~Fdone 
            Ftrind = max(find(rtime<=Ftrgs(s))) - ss(1) + 1;
            if (Ftrind+Fw) <= length(fltw)
                Fwin(:,s) = fltw(Ftrind:Ftrind+Fw)';
                if s == size(Ftrgs,2)
                    Fdone = 1;            
                    Fwin = Fwin(:,1:s-1);
                    Fsmu(:,j) = mean(Fwin,2);
                    Fcount(j) = s-1;   
                end    
            else
                Fdone = 1;
                Fwin = Fwin(:,1:s-1);
                Fsmu(:,j) = mean(Fwin,2);
                Fcount(j) = s-1;
            end
        end    
        
        if ~Pdone
            Ptrind = dsearchn(rtime',Ptrgs(s));
            if (Ptrind+Pw/2) <= length(emg)
                Pwin(:,s) = emg(j,round(Ptrind-Pw/2:Ptrind+Pw/2))';
                if s == size(Ptrgs,2)
                    Pdone = 1;            
                    Pwin = Pwin(:,1:s-1);
                    Psmu(:,j) = mean(Pwin,2);
                    Pcount(j) = s-1;   
                end    
            else
                Pdone = 1;            
                Pwin = Pwin(:,1:s-1);
                Psmu(:,j) = mean(Pwin,2);
                Pcount(j) = s-1;   
            end
        end    
        
        if Fdone & Pdone
            break           
        end    
    end        
    
end

SMUF = Fsmu;    SMUP = Psmu;

SMUFonset = find(ttime==0);
baseline = SMUF(SMUFonset,:);

bSMUF = SMUF - repmat(baseline,[size(SMUF,1) 1]);
mSMUF = mean(bSMUF,2);
mSMUP = mean(SMUP,2);

figure(12);	clf

subplot(223), cla
plot(ttime,SMUF);   hold on;
plot(ttime,mSMUF,'k','LineWidth',2);    axis tight;
set(gca, 'XLim', [ttime(1) 100]);
xlabel('Time (ms)');
ylabel('Twitches (arbitrary units)');

subplot(224), cla
plot(etime,SMUP);   hold on;
plot(etime,mSMUP,'k','LineWidth',2);    axis tight;
xlabel('Time (ms)');
ylabel('SMUPs (arbitrary units)');

for j = 1 : p.nsta
    
    disp(['Displaying results of STA number ' num2str(j) ' of ' num2str(p.nsta) '.']);

    tw = twitch(0, ttime, p.twtforce(p.sMUs(j)), p.tc(p.sMUs(j)));
    subplot(221), cla
    plot(ttime, tw, ':'); hold on;
    plot(ttime, bSMUF(:,j)');
    axis tight; title(['STA Force: Sample Unit ' num2str(j)]);
    xlabel('Time (ms)');
    ylabel('Twitch (arbitrary units)');

    sm = muap2(p.del, etime, p.twtforce(p.sMUs(j)));
    subplot(222), cla
    plot(etime, sm, ':'); hold on;
    plot(etime, SMUP(:,j)');
    axis tight; title(['STA EMG: Sample Unit ' num2str(j)]);
    xlabel('Time (ms)');
    ylabel('SMUP (arbitrary units)');
    
    if j ~= p.nsta
        disp('dbcont to proceed to next STA plot.')
        keyboard
    end    
            
end            