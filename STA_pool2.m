function [ltw,emg,sMU,frate,trigs,smuFs,smuPs] = STA_pool2(excite,p,rtime,numspk,ntrigs)
%%--------------------------------------------------------
%% this function simulates a motor neuron pool 
%% with Gaussian noise on the isi distribution
%% it also selects a sample motor unit for triggering STA
%% inputs are:
%%       excite = excitation level (scalar 0-1)
%%       p = params of pool
%%       rtime = time vector
%%       numspk = maximum number of spikes in a spike train
%%       ntrigs = maximum number of triggers (based on low rate requirement)
%% outputs are:
%%       ltw = muscle force time series
%%       emg = emg time series
%%       sMU = index of the selected sample motor unit
%%       frate = firing rate of the selected motor unit
%%       trigs = spike times of selected motor unit
%%       smuFs = force trace from selected MU alone
%%       smuPs = emg trace from selected MU alone

%% if OUT OF MEMORY, replace twitch.m with twitch_memory.m (line 73),
%% and replace muap2.m with muap2_memory.m (line 76)

drive = excite*p.maxe;
fr = (p.gain*(drive - p.rte) + p.mfr );
fr(drive<p.rte) = 0;  %% If drive, aka. total excitation, is less than threshold set firing to zero
fr(fr>p.pfr) = p.pfr(fr>p.pfr);  %% If firing rate is greater than peak rate set to peak  
gotit = 0;

while ~gotit
    candidates = find( (fr<p.frmax) & (fr>0) );
    if isempty(find(fr>0))
        error(['The excitation is too low - there are no active units.']);
    elseif isempty(find(fr<p.frmax))
        error(['The excitation is too high - there are no slow-firing units.']);
    end
    
    candies_shuffled = candidates(randperm(length(candidates)));
    sMU = candies_shuffled(1);

    if ~ismember(sMU,nonzeros(p.sMUs))
        gotit = 1;
    else
        disp(['MU # ' num2str(sMU) ' tried to be sampled more than once!']);
        if all(ismember(candidates,nonzeros(p.sMUs)));
            disp('All the candidate MUs have been sampled already - restarting this contraction at a different level.');
            sMU = 0;    frate = 0;  trigs = zeros(1,ntrigs);
            ltw = zeros(1,length(rtime));   emg = ltw;
            fltw = ltw;
            smuFs = ltw;    smuPs = emg;
            return
        end    
    end    
end

frate = fr(sMU);

twlin = zeros(length(rtime),p.n);   emglin = zeros(length(rtime),p.n);

%% sim for each MN
for i = 1 : p.n
    
    if (fr(i)>0)  %if this unit is active
        
        disp(['Calculating contribution from unit number ' num2str(i) '.']);
        misi = 1000/fr(i);                      %% mean inter-spike interval
        
        %% set ISI distribution
        sisi = misi*p.cv;                       %% get std of ISI dist
        disi = normrnd(misi,sisi,1,numspk);     %% generate an ISI dist
        disi = cumsum(disi);                    %% generate spike times
        
        %% do twitch function
        twlin(:,i) = twitch(disi,rtime,p.twtforce(i),p.tc(i));
        %% do muap function
        emglin(:,i) = muap2(disi,rtime,p.twtforce(i));

        %% calc non-linear force gain as in Fuglevand
        rate = p.tc(i)/misi;
        sc = (1-exp(-2*0.4.^3))/0.4;            % gain value predicted for a norm. rate of 0.4
        rg = (1-exp(-2*(rate.^3)))./(rate*sc);
        rg(rate<0.4) = 1;
        twlin(:,i) = twlin(:,i)*rg;             %% get matrix of all twitches with gain

        if i == sMU
            trigs = disi(1:ntrigs);
        end

    end
    
    pause(0.1)
    
end
            
if nargin > 5
    smuFs = twlin(:,sMU)';
    smuPs = emglin(:,sMU)';
end    

ltw = sum(twlin,2)';
emg = sum(emglin,2)';

return