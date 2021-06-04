function par=sobie_etal_2002_combined_parameters_3(cc0,csr0)

par.ATPtotal=455; % uM
par.konATP=0.225; % uM-1ms-1
par.koffATP=45; % ms-1
par.CaMfree=24; % uM
par.koffCaM=0.238; % ms-1
par.konCaM=0.034; % uM-1ms-1
par.CaMtotal=24; %uM

% par.Vcell=4.4e-19;%m^3
% par.Vcell=8.5770e-21;
par.Vcell=8.5770e-3/6;%um^2
par.konTnC=0.0327;
par.koffTnC=0.0196;
par.TnCtotal=70;
par.konF4=0.1;
par.koffF4=0.11;
par.F4total=24;

par.konSRbuf=0.115;
par.koffSRbuf=0.1;
par.SRbuftotal=0;
par.numserca=1;
par.Vserca=200;
par.Kserca=0.184;
par.mserca=4;

par.tau=csr0*(par.Kserca^par.mserca+cc0^par.mserca)/(par.numserca*par.Vserca*cc0^par.mserca);

% Walker's (cyted Tran et al) sercas model parameters:
par.Kdi=910; %uM
par.Kdsr=2240; %uM
par.Ap=150; %uM

Kij=(cc0/par.Kdi).^2;
Ksr=(csr0/par.Kdsr).^2;
vcycle=(3.24873*(10^12)*(Kij.^2)+Kij.*(9.17846*(10^6)-11478.2.*Ksr)-0.329904.*Ksr)./...
    (0.104217+17.293.*Ksr+Kij.*(1.75583.*(10^6)+7.61673.*(10^6).*Ksr)+(Kij.^2).*(10^11).*(6.08462+4.50544.*Ksr));

par.tau_walker=csr0./(2*vcycle*par.Ap);


par.F=9.6485e4;
par.Vdyad=1.71e-21;
par.Vdyad2=1.71e-3/0.2;%um^2
% par.Vdyad2=1;
par.gcyto=160;
par.gryr=0.033;
par.kcoop=1;
par.gillc=2e-16;
par.tlcc=0.5;
par.lag_lcc=2; %ms
par.t_lcc=0.5; %ms
% par.lag_lcc=0;
% par.t_lcc=5;

par.koffSLbuf=1;
par.konSLbuf=0.115;
par.SLbuftotal=0;


% % par.grefill=0.04;
par.grefill=0.095;
par.CSQtotal=10000;
par.KCSQ=800;

par.Dcyto=0.22;
par.Ddyad=0.22;
par.Dnsr=0.0733;
par.DF4Ca=0.042;
par.DF4=0.042;
par.DCaM=0.025;
par.DCaCaM=0.025;
par.DATP=0.14;
par.DCaATP=0.14;

par.numdiad=1;