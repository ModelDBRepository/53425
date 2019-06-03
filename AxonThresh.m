function Is = AxonThresh(depth,diam)

% calculate axon thresholds based on all combinations of depth and diameter
% depth: length M
% diam: length N
% Is: size MxN

L = 400;                                % length of axon field (cm)
dVmt = 10;                              % depolarization threshold (mV)
dAC = 2.5;                              % distance between anode and cathode (cm)
dz = 0.01;                              % axial spatial resolution (cm)

zA = L/2 + dAC/2;                       % anode position (cm)
zC = zA - dAC;                          % cathode position (cm)
z = 0:dz:L;                             % vector of axial positions (cm)

Rm = 5000;                    % specific membrane resistance (ohm.cm^2)
Ra = 70;                      % specific axial resistance (ohm.cm)
re = 350;                     % extracellular resistance (ohm.cm)

Is = zeros(length(depth),length(diam));

for j = 1: length(diam)
    
    r = diam(j) / 1e4 / 2;                  % axon radius (cm)
    rm = Rm./(2*pi.*r);                     % membrane resistance (ohm.cm)
    ra = Ra./(pi.*r.^2);                    % axial resistance (ohm/cm)
    k = sqrt(ra/rm);                        % spatial frequency of intracellular potential
                                            % 1/um
    for i = 1 : length(depth)
        
        d = depth(i) / 1e4;                     % convert depth to cm
    
        dA = sqrt( (d-r)^2 + (zA-z).^2 );       % absolute distances from electrodes
        dC = sqrt( (d-r)^2 + (zC-z).^2 );       % to axon surface coordinates (cm)

        VA = re./(4*pi*dA);                     % potential outside membrane due to 1mA anode (mV)
        VC = re./(4*pi*dC);                     % potential outside membrane due to 1mA cathode (mV)
        Ve = VA - VC;                           % potential applied to membrane (mV)

        dVi = Visolve(Ve,rm,ra,dz);
        dVm1 = max(dVi - Ve);

        Is(i,j) = dVmt ./ dVm1;
        
    end
    
end    

return