function [m42_new,h42_new,m24_new,h24_new]=update_ip3r_states(cm,m42,h42,m24,h24,lambda_h42,ip3r_par,par)
m42_inf=cm.^ip3r_par.n42./(ip3r_par.k42.^ip3r_par.n42+cm.^ip3r_par.n42);

m42_new=m42+par.dt*ip3r_par.lambda_m42.*(m42_inf-m42);
% YK1=par.dt*ip3r_par.lambda_m42.*(m42_inf-m42);
% YK2=par.dt*ip3r_par.lambda_m42.*(m42_inf-(m42+YK1./2));
% YK3=par.dt*ip3r_par.lambda_m42.*(m42_inf-(m42+YK2./2));
% YK4=par.dt*ip3r_par.lambda_m42.*(m42_inf-(m42+YK3));
% m42_new=m42+(YK1+2*YK2+2*YK3+YK4)/6;

h42_inf=ip3r_par.kn42.^ip3r_par.nn42./(ip3r_par.kn42^ip3r_par.nn42+cm.^ip3r_par.nn42);

h42_new=h42+par.dt*lambda_h42.*(h42_inf-h42);
% YK1=par.dt*lambda_h42.*(h42_inf-h42);
% YK2=par.dt*lambda_h42.*(h42_inf-(h42+YK1/2));
% YK3=par.dt*lambda_h42.*(h42_inf-(h42+YK2/2));
% YK4=par.dt*lambda_h42.*(h42_inf-(h42+YK3));
% h42_new=h42+(YK1+2*YK2+2*YK3+YK4)./6;

m24_inf=cm.^ip3r_par.n24./(ip3r_par.k24^ip3r_par.n42+cm.^ip3r_par.n24);

m24_new=m24+par.dt*ip3r_par.lambda_m24.*(m24_inf-m24);
% YK1=par.dt*ip3r_par.lambda_m24.*(m24_inf-m24);
% YK2=par.dt*ip3r_par.lambda_m24.*(m24_inf-(m24+YK1./2));
% YK3=par.dt*ip3r_par.lambda_m24.*(m24_inf-(m24+YK2./2));
% YK4=par.dt*ip3r_par.lambda_m24.*(m24_inf-(m24+YK3));
% m24_new=m24+(YK1+2*YK2+2*YK3+YK4)./6;

h24_inf=ip3r_par.kn24.^ip3r_par.nn24./(ip3r_par.kn24.^ip3r_par.nn24+cm.^ip3r_par.nn24);

h24_new=h24+par.dt*ip3r_par.lambda_h24.*(h24_inf-h24);
% YK1=par.dt*ip3r_par.lambda_h24.*(h24_inf-h24);
% YK2=par.dt*ip3r_par.lambda_h24.*(h24_inf-(h24+YK1/2));
% YK3=par.dt*ip3r_par.lambda_h24.*(h24_inf-(h24+YK2/2));
% YK4=par.dt*ip3r_par.lambda_h24.*(h24_inf-(h24+YK3));
% h24_new=h24+(YK1+2*YK2+2*YK3+YK4)./6;

end