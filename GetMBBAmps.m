% script to keep track of Revised Incremental Stimulation algorithm parameters
n = stepcnt;

AFall = max(ltw,[],2);          [AFu, AFui] = unique(AFall);
data.AFmx(n) = max(AFall);      data.AFmn(n) = min(AFall);
APall = max(emg,[],2);          [APu, APui] = unique(APall);
data.APmx(n) = max(APall);      data.APmn(n) = min(APall);

if n == 1 
    nzswi = min(find(AFall));
    AFall = nonzeros(AFall);    AFu = nonzeros(AFu);    AFui = AFui(find(nrec(AFui)));
    APall = nonzeros(APall);    APu = nonzeros(APu);    APui = APui(find(nrec(APui)));
else
	nzswi = 1;
end   

if length(AFu) == 1
    data.AFpd(n) = AFu;   
    data.Fswi(n) = nzswi;
    FPDcount = length(AFall);
    data.APpd(n) = APu;
    data.Pswi(n) = nzswi;
    PPDcount = length(APall);
    if FPDcount < ceil(p.sweeps/2)
        data.AFpd(n) = 0;
        data.APpd(n) = 0;
        %%% ERROR CATCHING %%%
            if PPDcount >= ceil(p.sweeps/2)
                error('There were not enough force responses, but there were enough emg???');
            end    
        %%% ERROR CATCHING %%%
    end    
    %%% ERROR CATCHING %%%
        if length(APu) ~= 1
            error('There was only one force response level, but more than one emg response level.');
        end    
    %%% ERROR CATCHING %%%
else
    
    FPDcount = hist(AFall,AFu);
    FPDi = find(FPDcount == max(FPDcount));
    data.AFpd(n) = min(AFu(FPDi));
    data.Fswi(n) = AFui(min(FPDi));
    
    PPDcount = hist(APall,APu);
    PPDi = find(PPDcount == max(PPDcount));
    data.APpd(n) = min(APu(PPDi));
    data.Pswi(n) = APui(min(PPDi));

end