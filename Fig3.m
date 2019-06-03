%% ModelDB 

%%--------------------------------------------------------
%% J. Neural Eng. 2 (2005) 17–34
%% "Simulations of motor unit number estimation techniques"
%% AUTHOR: Lora A Major and Kelvin E Jones
%% SUBMITTED:  19NOV2004
%% FIGURE 3
%% MODIFIED DATE: 30MAY2005
%%--------------------------------------------------------

%% define axon diameters (um)
Diams = 4 : 0.5: 16;

%% define axon depths (um) from fascicle survace + 2 mm to skin surface
Depths = [200 : 10 : 400] + 2000;

%% calculate axon tresholds
Is = AxonThresh(Depths,Diams);

figure(3); clf

Is = AxonThresh(Depths, Diams);
surf(Depths, Diams, Is');   axis tight;     axis ij
set(gca, 'YTick', 5:5:15);
view(50, 25)
xlabel('Axon Depth (um)');
ylabel('Axon Diameter (um)');
zlabel('Axon Threshold (mA)');