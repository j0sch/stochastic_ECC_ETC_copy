function Yn_plus_1=RK4_Sobie_Cao_2_diffusion(Yn,par,ip3r_par,RyR_open,ip3r_open,simt)
% Solves Ca signalling ODEs using 4th order Runge-Kutta method.
% Matrix Yn contains system variables at time tn
% par - structure containing all parameters of the model
% RyR_open - matrix containing information about number of open RyR
% receptors at each discretization point.

% Function Sobie_fluxes.m contains discretized fluxes with known states of
% each receptor
K1=Sobie_Cao_fluxes_combined_2(Yn,par,ip3r_par,RyR_open,ip3r_open,simt);
Yn_plus_1=Yn+par.dt*K1;
% K2=Sobie_Cao_fluxes_combined_2(Yn+par.dt*K1/2,par,ip3r_par,RyR_open,ip3r_open,simt);
% K3=Sobie_Cao_fluxes_combined_2(Yn+par.dt*K2/2,par,ip3r_par,RyR_open,ip3r_open,simt);
% K4=Sobie_Cao_fluxes_combined_2(Yn+par.dt*K3,par,ip3r_par,RyR_open,ip3r_open,simt);
% 
% Yn_plus_1=Yn+(par.dt/6)*(K1+2*K2+2*K3+K4);

end
