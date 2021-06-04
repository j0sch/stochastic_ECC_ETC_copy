function Y=Sobie_Cao_fluxes_combined_2(Y0,par,par_ip3r,RyR_open,IP3_open,simt)
% Function calculating main fluxes
% Y0=[Cai, Cajsr, Cansr, CaM, CaCaM, ATP, CaATP, F4, F4Ca, TnC];
%      1     2      3     4     5     6     7     8    9    10   

vs=size(Y0,1);

JCaM=(par.koffCaM*Y0(:,5)-par.konCaM*Y0(:,1).*Y0(:,4));
JATP=(par.koffATP*Y0(:,7)-par.konATP*Y0(:,1).*Y0(:,6));
JF4= (par.koffF4*Y0(:,9)-par.konF4*Y0(:,1).*Y0(:,8));
Jtnc=(par.koffTnC*(par.TnCtotal-Y0(:,10))-par.konTnC*Y0(:,1).*Y0(:,10)).*(~par.positions);

% Diffusion terms: 
J_diff=zeros(vs,10);
% 
J_diff(2:end-1,:)=(Y0(3:end,:)+Y0(1:end-2,:)-2*Y0(2:end-1,:))./(par.dx^2);
J_diff(1,:)=2.*(Y0(2,:)-Y0(1,:))./(par.dx^2);
J_diff(end,:)=2.*(Y0(end-1,:)-Y0(end,:))./(par.dx^2);


% % 
% % Walker's (cyted Tran et al) sercas:
Kij=(Y0(:,1)./par.Kdi).^2;
Ksr=(Y0(:,3)./par.Kdsr).^2;
vcycle=(3.24873*(10^12)*(Kij.^2)+Kij.*(9.17846*(10^6)-11478.2.*Ksr)-0.329904.*Ksr)./...
    (0.104217+17.293.*Ksr+Kij.*(1.75583.*(10^6)+7.61673.*(10^6).*Ksr)+(Kij.^2).*(10^11).*(6.08462+4.50544.*Ksr));
Jserca=par.scale_sercas*(2*vcycle*par.Ap-Y0(:,3)./par.tau_walker).*par.sercas;%(~par.positions);
% 
% Jserca=(2*vcycle*par.Ap).*par.sercas;%(~par.positions);


% Jserca=Jserca.*(Jserca>0);
% 
% %  VJ model +leakage
% Jserca=((((par.numserca*(par.Vserca*Y0(:,1).^par.mserca))./...
%     (par.Kserca.^par.mserca+Y0(:,1).^par.mserca)))-Y0(:,3)/par.tau).*(~par.positions);
% % % % % No sercas
% % % % Jserca=0;z

Jlcc=0;
% Jlcc=(par.gillc*1000./(2*par.F*par.Vdyad*1)).*(simt>par.lag_lcc).*(simt<=(par.lag_lcc+par.t_lcc));
Jrel=(par.gryr*RyR_open+par_ip3r.kipr*IP3_open).*(Y0(:,2)-Y0(:,1));
% Jrel=(par.gryr*RyR_open).*(Y0(:,2)-Y0(:,1));
% JRyR=0;

beta_jsr=(1+(par.CSQtotal*par.KCSQ)./((Y0(:,2)+par.KCSQ).^2)).^-1;
% beta_jsr=1;
Jrefill=(Y0(:,3)-Y0(:,2)).*par.grefill.*par.positions;
% Jrefill=(Y0(:,2)-Y0(:,1)).*par.grefill;



% Y=[par.Dcyto*J_diff(:,1)+JF4+JCaM+JATP+Jtnc-Jserca+(Jrel+Jlcc).*(par.positions),...
%     beta_jsr.*(Jrefill-(par.Vdyad*1e18*5/par.Vjsr).*Jrel),...
%     par.Dnsr*J_diff(:,3)+(Jserca*((par.Vcell)/par.Vnsr)-Jrefill*par.Vjsr/par.Vnsr),...% 0.*(Jserca-Jrefill),...%  
%     JCaM+par.DCaM*J_diff(:,4),...
%     -JCaM+par.DCaCaM*J_diff(:,5),...
%     JATP+par.DATP*J_diff(:,6),...
%     -JATP+par.DCaATP*J_diff(:,7),...
%     JF4+par.DF4*J_diff(:,8),...
%     -JF4+par.DF4Ca*J_diff(:,9),...
%     Jtnc];


Y=[par.Dcyto*J_diff(:,1)+JF4+JCaM+JATP+Jtnc-Jserca+(Jrel+Jlcc).*(par.positions),...
    beta_jsr.*(Jrefill-Jrel),... %  par.Dnsr*J_diff(:,2)+ for now to test possibility. Generally region should be smaller, thus boundaries are not the same as this case
    par.Dnsr*J_diff(:,3)+(Jserca-Jrefill),...% 0.*(Jserca-Jrefill),...%  
    JCaM+par.DCaM*J_diff(:,4),...
    -JCaM+par.DCaCaM*J_diff(:,5),...
    JATP+par.DATP*J_diff(:,6),...
    -JATP+par.DCaATP*J_diff(:,7),...
    JF4+par.DF4*J_diff(:,8),...
    -JF4+par.DF4Ca*J_diff(:,9),...
    Jtnc];
