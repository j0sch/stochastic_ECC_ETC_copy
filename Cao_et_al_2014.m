%% Written by Agne Tilunaite 2019
%% Some comments added by Hilary Hunt 2020
% This function outputs the parameter values to the Cao et al 2014 IP3R
% model.

% The meaning of each parameter is described in the 2014 PLoS Comp Bio
% paper.
% To increase the upper range of [Ca2+] that will produce IP3R activity,
% increasing kn42 and/or decreasing kn24 may help.
function par_ip3r=Cao_et_al_2014(Num_ipr,c0,ip3c)

% par_ip3r.p=0.15;
% %  Change IP3 concentration
par_ip3r.p=ip3c;
par_ip3r.v42=110*par_ip3r.p^2./(par_ip3r.p^2+0.1^2);%
par_ip3r.k42=0.49+0.543*par_ip3r.p^3./(par_ip3r.p^3+4^3);%
par_ip3r.n42=3;%
par_ip3r.kn42=0.41+25*par_ip3r.p.^3./(par_ip3r.p.^3+6.5^3);%
par_ip3r.nn42=3;%
par_ip3r.a42=1.8*par_ip3r.p.^2./(par_ip3r.p^2+0.58^2);%
par_ip3r.v24=62+880./(par_ip3r.p^2+4);%
par_ip3r.k24=0.35;%
par_ip3r.n24=3;%
par_ip3r.kn24=80;%
par_ip3r.nn24=2;%
par_ip3r.a24=1+5./(par_ip3r.p^2+0.5^2);%

par_ip3r.scale=(4010/(4010+10500));


par_ip3r.kipr=0.05/(Num_ipr*(1-par_ip3r.scale));

par_ip3r.h24_inf=par_ip3r.kn24^par_ip3r.nn24./...
    (par_ip3r.kn24^par_ip3r.nn24+c0^par_ip3r.nn24);%

par_ip3r.lambda_h24=40;%

par_ip3r.m24_inf=c0^par_ip3r.n24./(par_ip3r.k24^par_ip3r.n24+c0^par_ip3r.n24);%

par_ip3r.lambda_m24=100;%

par_ip3r.m42_inf=c0^par_ip3r.n42./(par_ip3r.k42^par_ip3r.n42+c0^par_ip3r.n42);%

par_ip3r.lambda_m42=100;%

par_ip3r.h42_inf=par_ip3r.kn42^par_ip3r.nn42/(par_ip3r.kn42^par_ip3r.nn42+c0^par_ip3r.nn42);%

par_ip3r.ah42=0.5;
par_ip3r.vh42=20;

