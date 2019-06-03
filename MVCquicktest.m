disp('Calculating a table of force output for 0 to 100% excitation.');

nonlin = 1;             %% include force-rate nonlinearity (0 of not)

%% define time scale
p.dt = 0.5;              %% Time step (0.5 ms -> 2000kHz sampling rate)
ftime = 0:p.dt:4000;

p.lev = 0.02:0.01:1.0;        %% levels of excitation to test
nlev = length(p.lev);

mF = zeros(nlev,1);

for j = 1 : nlev
        drive = p.lev(j)*p.maxe;
        fr = (p.gain*(drive - p.rte) + p.mfr );
        fr(drive<p.rte) = 0;  %% If drive, aka. total excitation, is less than threshold set firing to zero
        fr(fr>p.pfr) = p.pfr(fr>p.pfr);  %% If firing rate is greater than peak rate set to peak  

        F = zeros(p.n,1);

        %% sim for each MN
        for i = 1 : p.n
    
            if (fr(i)>0)  %if this unit is active
        
                %% integrate twitch function
                F(i) = fr(i)/1000 * sum(twitch(0, ftime, p.twtforce(i), p.tc(i))) * p.dt;       
            
                if nonlin
                    %% calc non-linear force gain as in Fuglevand
                    rate = p.tc(i)*fr(i)/1000;
                    sc = (1-exp(-2*0.4.^3))/0.4;            % gain value predicted for a norm. rate of 0.4
                    rg = (1-exp(-2*(rate.^3)))./(rate*sc);
                    rg(rate<0.4) = 1;
                    F(i) = F(i) * rg;
                end   
        
            end
    
        end
            
        mF(j) = sum(F);
        
end

lev = p.lev;

clear nonlin nlev drive fr F rate sc rg