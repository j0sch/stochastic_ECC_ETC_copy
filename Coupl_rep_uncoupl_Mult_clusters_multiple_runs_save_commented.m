%% Originally Coupl_rep_uncoupl_Mult_clusters_multiple_runs_save.m
%% Written by Agne Tilunaite 2019
%% Edited by Hilary Hunt 2020
numSimulations=100;
c0=0.1;
% Length of simulation
simLength=3000; %ms
% Distance between clusters
clustDist=1.8; %um
% Frequency of triggered sparks
triggerFreq=1000;

casr0=800*ones(numSimulations,1);
% Change scale to 10 for IP3R type II
scale=1;
for numIP3R=[0,5,10,20]
% Define parameters
% Use parfor
parfor fop=1:numel(casr0)
             
              
    % Get parameters from Sobie 2002   
    % first param specifies initial cytosolic calcium (uM)
    % second param specifies initial SR calcium (uM)
    par=sobie_etal_2002_combined_parameters_3(c0,casr0(fop));
    par.casr_0=casr0(fop);
    par.sc=10;
    par.scale_sercas=1;
    par.gillc=0;
    par.gryr=0.02;
    par.dt0=1e-4; %ms
    par.dx=0.2/5; %um
  
    par.dist=clustDist; % um distance between coupled and uncoupled clusters
    par.num_of_fixed_ip3r=10;
    par.num_of_varied_ip3r=numIP3R;
    
    par.ryr_op_times=triggerFreq;
    par.keep_ryr_open=5;
   
    % Get Cao 2014 parameters 
    % first param specifies number of IP3Rs
    % second param specifies initial cytosolic calcium concentration at
    % IP3Rs
    % third param specifiec ip3 conc
    ip3r_par=Cao_et_al_2014(1,c0*scale,0.15);
    ip3r_par.kipr=0.02;
     ip3r_par.rescale_ca=scale;
%     par.kip3r=ip3r_par.kipr;

    linescan_length=ceil((2+par.dist+par.dist+2)/par.dx);

    FOV=[1,linescan_length];

    RyR_place=ones(FOV(1),FOV(2)); % In other words RyR_place=RyR_matrix>0

    % Storage of each receptor (RyR and then IP3R) behaviour and its transition times:
    % Flatten information about all receptors:
    % 1st column - receptor position in the flattened receptor matrix
    % 2nd column - the state of the receptor
    % 3rd column - waiting time till receptor changes state
    % 4th column - random number to determine receptor transition times

    numofryrs=90*3;
    numofip3rs=par.num_of_fixed_ip3r*2+2*par.num_of_varied_ip3r;

    total_rec=numofryrs+numofip3rs;


    RyR_release_place=ceil(linescan_length/2); % for uncoupled one
    RyR_sparse=zeros(numofryrs+numofip3rs,4);
    % RyR_pos=RyR_release_place*ones(numofryrs,1);

    RyR_open_sums=zeros(linescan_length,1);
    IP3R_open_sums=zeros(linescan_length,1);

    initiationT=simLength;
    initiationiter=ceil(initiationT/par.dt0);
    titer=ceil(simLength/par.dt0)/100;

    % indx_fl=0;
    % for aa=1:numel(RyR_pos)
    %    RyR_sparse(indx_fl+1:indx_fl+1,1)=RyR_pos(aa);
    %    indx_fl=indx_fl+1;
    % end

%     RyR_sparse(1,1)=RyR_release_place-(par.dist/par.dx)-1;       %1st coupled cluster
%     RyR_sparse(2,1)=RyR_release_place-(par.dist/par.dx)+1;

% %     RyR_sparse(1,1)=1;       %1st coupled cluster
% %     RyR_sparse(2,1)=3;
% %     
% %     RyR_sparse(3:45+2)=RyR_release_place-1; % uncoupled cluster
% %     RyR_sparse(45+3:90+2)=RyR_release_place+1;
% %     
% %     RyR_sparse(93)=RyR_release_place+(par.dist/par.dx)-1; % 2nd coupled cluster
% %     RyR_sparse(94)=RyR_release_place+(par.dist/par.dx)+1;
    
    div=numofryrs/6;

    RyR_sparse(1:div,1)=RyR_release_place-(par.dist/par.dx)-1;       %1st coupled cluster
    RyR_sparse(div+1:2*div,1)=RyR_release_place-(par.dist/par.dx)+1;
    
    RyR_sparse(2*div+1:3*div)=RyR_release_place-1; % uncoupled cluster
    RyR_sparse(3*div+1:4*div)=RyR_release_place+1;
    
    RyR_sparse(4*div+1:5*div)=RyR_release_place+(par.dist/par.dx)-1; % 2nd coupled cluster
    RyR_sparse(5*div+1:numofryrs)=RyR_release_place+(par.dist/par.dx)+1;
    
    
    RyR_sparse(1+numofryrs:1+numofryrs+par.num_of_fixed_ip3r,1)=RyR_release_place-(par.dist/par.dx);
% %     RyR_sparse(1+numofryrs:1+numofryrs+par.num_of_fixed_ip3r,1)=2;
    RyR_sparse(1+numofryrs+par.num_of_fixed_ip3r:numofryrs+par.num_of_fixed_ip3r+par.num_of_varied_ip3r,1)=RyR_release_place-2;
    RyR_sparse(1+numofryrs+par.num_of_fixed_ip3r+par.num_of_varied_ip3r:numofryrs+numofip3rs-par.num_of_fixed_ip3r,1)=RyR_release_place+2;
    RyR_sparse(1+numofryrs+numofip3rs-par.num_of_fixed_ip3r:total_rec,1)=RyR_release_place+(par.dist/par.dx);

    par.positions=zeros(linescan_length,1);
    par.positions(RyR_release_place-2:RyR_release_place+2)=1;
% %     par.positions(1:3)=1;
    par.positions(RyR_release_place-(par.dist/par.dx)-2:RyR_release_place-(par.dist/par.dx)+2)=1;
    par.positions(RyR_release_place+(par.dist/par.dx)-2:RyR_release_place+(par.dist/par.dx)+2)=1;

    par.sercas=~par.positions;
     
    %State 1 closed; state 2 open
    Ch_mx=[-1;1];

    RyR_sparse(:,2)=1; % state of each receptor
% %     RyR_sparse(3:33,2)=2; % state of each receptor
% %     RyR_sparse(48:48+30,2)=2; % state of each receptor

% Change the following two lines to '=1' instead of '=2' to start the
% simulation without triggered sparks at the two coupled dyads.
    RyR_sparse(1:30,2)=2; % state of each receptor
    RyR_sparse(div+1:div+30,2)=2; % state of each receptor
    RyR_sparse(4*div+1:4*div+30,2)=2; % state of each receptor
    RyR_sparse(5*div+1:5*div+30,2)=2; % state of each receptor
    
    
    RyR_open_sums(RyR_sparse(1,1))=sum(RyR_sparse(1:div,2)==2);
    RyR_open_sums(RyR_sparse(div+1,1))=sum(RyR_sparse(div+1:2*div,2)==2);
    RyR_open_sums(RyR_sparse(2*div+1,1))=sum(RyR_sparse(2*div+1:3*div,2)==2);
    RyR_open_sums(RyR_sparse(3*div+1,1))=sum(RyR_sparse(3*div+1:4*div,2)==2);
    RyR_open_sums(RyR_sparse(4*div+1,1))=sum(RyR_sparse(4*div+1:5*div,2)==2);
    RyR_open_sums(RyR_sparse(5*div+1,1))=sum(RyR_sparse(5*div+1:numofryrs,2)==2);

% %     RyR_open_sums(RyR_sparse(1,1))=RyR_sparse(1,2)==2;
% %     RyR_open_sums(RyR_sparse(2,1))=RyR_sparse(2,2)==2;
% %     RyR_open_sums(RyR_sparse(3,1))=sum(RyR_sparse(3:45+2,2)==2);
% %     RyR_open_sums(RyR_sparse(48,1))=sum(RyR_sparse(45+3:90+2,2)==2);
% %     RyR_open_sums(RyR_sparse(93,1))=RyR_sparse(93,2)==2;
% %     RyR_open_sums(RyR_sparse(94,1))=RyR_sparse(94,2)==2;

    RyR_sparse(:,4)=rand(total_rec,1);

    % Variables I will need to track

    % Y0=[Cai, Cajsr, Cansr, CaM, CaCaM, ATP, CaATP, F4, F4Ca, TnC];
    %      1     2      3     4     5     6     7     8    9    10   

    % RyR &SR related regions

   

    RyR_open_save=zeros(linescan_length,titer);
    IP3R_open_save=zeros(linescan_length,titer);
    Cai_save=zeros(linescan_length,titer);
    Cajsr_save=zeros(linescan_length,titer);
    Caf4_save=zeros(linescan_length,titer);
    ATP_save=zeros(linescan_length,titer);
    CaM_save=zeros(linescan_length,titer);
    TpC_save=zeros(linescan_length,titer);
    Cansr_save=zeros(linescan_length,titer);
    beta_check=zeros(linescan_length,titer);
    time=zeros(1,titer);

    % initial conditions:
    Cai=c0.*ones(linescan_length,1);

    Cajsr=par.casr_0*ones(linescan_length,1);
    Cansr=par.casr_0*ones(linescan_length,1);
    CaM0=par.koffCaM.*par.CaMtotal/(par.koffCaM+par.konCaM*0.1);
    ATP0=par.koffATP.*par.ATPtotal/(par.koffATP+par.konATP*0.1);
    F40=par.koffF4.*par.F4total/(par.koffF4+par.konF4*0.1);
    TnC0=par.koffTnC*par.TnCtotal/(par.koffTnC+par.konTnC*0.1);

    CaM=CaM0*ones(linescan_length,1);
    CaCaM=(par.CaMtotal-CaM0)*ones(linescan_length,1);
    ATP=ATP0*ones(linescan_length,1);
    CaATP=(par.ATPtotal-ATP0)*ones(linescan_length,1);
    CaF4=(par.F4total-F40)*ones(linescan_length,1);
    F4=F40*ones(linescan_length,1);
    TnC=TnC0*ones(linescan_length,1);

    h24_inf=ip3r_par.kn24^ip3r_par.nn24./(ip3r_par.kn24^ip3r_par.nn24+(c0*scale).^ip3r_par.nn24);
    h24=h24_inf*ones(numofip3rs,1);

    m24_inf=c0.^ip3r_par.n24./(ip3r_par.k24^ip3r_par.n24+(c0*scale).^ip3r_par.n24);
    m24=m24_inf*ones(numofip3rs,1);

    m42_inf=c0.^ip3r_par.n42./(ip3r_par.k42.^ip3r_par.n42+(c0*scale).^ip3r_par.n42);
    m42=m42_inf*ones(numofip3rs,1);

    h42_inf=ip3r_par.kn42^ip3r_par.nn42./(ip3r_par.kn42^ip3r_par.nn42+(c0*scale).^ip3r_par.nn42);
    h42=h42_inf*ones(numofip3rs,1);

    tic
    par.dt=par.dt0;
    simt=0;
    rep=1;
    old_indx=0;
    for iter=1:titer
        dumt=0;
        while dumt<=(par.dt0*100)
           
            if (simt+dumt)>(rep*par.ryr_op_times) && (simt+dumt)<(rep*par.ryr_op_times+par.keep_ryr_open)
                % Trigger sparks at coupled dyads
                RyR_sparse(1:30,2)=2; % state of each receptor
                RyR_sparse(div+1:div+30,2)=2; % state of each receptor
                RyR_sparse(4*div+1:4*div+30,2)=2; % state of each receptor
                RyR_sparse(5*div+1:5*div+30,2)=2; % state of each receptor


                RyR_open_sums(RyR_sparse(1,1))=sum(RyR_sparse(1:div,2)==2);
                RyR_open_sums(RyR_sparse(div+1,1))=sum(RyR_sparse(div+1:2*div,2)==2);
                RyR_open_sums(RyR_sparse(2*div+1,1))=sum(RyR_sparse(2*div+1:3*div,2)==2);
                RyR_open_sums(RyR_sparse(3*div+1,1))=sum(RyR_sparse(3*div+1:4*div,2)==2);
                RyR_open_sums(RyR_sparse(4*div+1,1))=sum(RyR_sparse(4*div+1:5*div,2)==2);
                RyR_open_sums(RyR_sparse(5*div+1,1))=sum(RyR_sparse(5*div+1:numofryrs,2)==2);
                
                rep=ceil((simt+dumt)/par.ryr_op_times);

            end
            
           Yn=[Cai, Cajsr, Cansr, CaM,CaCaM, ATP, CaATP, F4, CaF4, TnC];
           Yn1=RK4_Sobie_Cao_2_diffusion(Yn,par,ip3r_par,RyR_open_sums,IP3R_open_sums,simt);             % values at the time step n+1


            % Dumm0 variable to find out real minimal time step between the state
            % change of each receptor.
% %             dt1=par.dt*ones(total_rec,1);
            dt1=par.dt0*ones(total_rec,1);



            %vectorised code again
            ryrVector = (1:numofryrs);
            ip3rVector=(numofryrs+1:total_rec);
            tanVector = [[zeros(numofryrs,1)...
                min(3.17*1e2*((Yn(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                max(0.250*((Yn(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)];...
                [zeros(numofip3rs,1)...
                ip3r_par.a42+ip3r_par.v42.*m42.*h42...
                (ip3r_par.a24+ip3r_par.v24.*(1-m24.*h24)).*ip3r_par.scale...
                zeros(numofip3rs,1)]];

            lambda_h42=(ip3r_par.vh42-ip3r_par.ah42)*(RyR_sparse(1+numofryrs:end,2)==2)+ip3r_par.ah42;

            cm=(120*(Yn1(RyR_sparse(ip3rVector),2)/200)-Yn1(RyR_sparse(ip3rVector),1)/10).*(RyR_sparse(1+numofryrs:end,2)==2)+Yn1(RyR_sparse(ip3rVector),1)/10;
            [m42_new,h42_new,m24_new,h24_new]=update_ip3r_states((cm*scale),m42,h42,m24,h24,lambda_h42,ip3r_par,par);

            tan1Vector = [[zeros(numofryrs,1)...
                min(3.17*1e2*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                max(0.250*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)];...
                [zeros(numofip3rs,1)...
                ip3r_par.a42+ip3r_par.v42.*m42_new.*h42_new...
                (ip3r_par.a24+ip3r_par.v24.*(1-m24_new.*h24_new)).*ip3r_par.scale...
                zeros(numofip3rs,1)]];

            tanIndices = [(RyR_sparse(:,2)*2)-1 RyR_sparse(:,2)*2];   
            tanVectorState = tanVector(([tanIndices-1]*total_rec)+[1:total_rec]');
            tan1VectorState = tan1Vector(([tanIndices-1]*total_rec)+[1:total_rec]');

            gOldVector = RyR_sparse(1:total_rec,3);

            gNewVector = gOldVector +...
                par.dt.*...
                (sum(tanVectorState,2)+sum(tan1VectorState,2))/2;

            xiVector = log(RyR_sparse(1:total_rec,4).^-1);  %% this is the best way to do element-wise inverse

            gNewGreaterThanXiBoolean = gNewVector > xiVector;

            dt1Vector = dt1; % I'm just defining this so I don't mess with your existing variables
            dt1Vector(gNewGreaterThanXiBoolean) = (xiVector(gNewGreaterThanXiBoolean) - gOldVector(gNewGreaterThanXiBoolean))./...
                (gNewVector(gNewGreaterThanXiBoolean) - gOldVector(gNewGreaterThanXiBoolean))*par.dt;

            [dt_min, indx_dt_min]=min(abs(dt1Vector));

%             disp([dt_min*1000,indx_dt_min,par.dt*1000])
    % % 
    % % %         dbstop if dt_min<=0
            if (dt_min==par.dt0) && (indx_dt_min==old_indx)
                indx_dt_min=0;
%                 disp([dt_min*1000,indx_dt_min,par.dt*1000])
            end
            old_indx=indx_dt_min;
                par.dt=dt_min;
    %             par.dt=max(dt_min,1e-4);

                % Recalculate system with minimal dt
                Yn1=RK4_Sobie_Cao_2_diffusion(Yn,par,ip3r_par,RyR_open_sums,IP3R_open_sums,simt);          % values at the time step n+1      
                if indx_dt_min==0
                    par.dt=par.dt0;
                    
                else
                   
%                     disp([dt_min*1000,indx_dt_min,par.dt*1000])
                    
    %             %replace the loop above with vectorised code
                ryrVector = (1:numofryrs);
                ip3rVector=(numofryrs+1:total_rec);

                tanVector = [[zeros(numofryrs,1)...
                    min(3.17*1e2*((Yn(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                    max(0.250*((Yn(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)];...
                    [zeros(numofip3rs,1)...
                    ip3r_par.a42+ip3r_par.v42.*m42.*h42...
                    (ip3r_par.a24+ip3r_par.v24.*(1-m24.*h24)).*ip3r_par.scale...
                    zeros(numofip3rs,1)]];


                lambda_h42=(ip3r_par.vh42-ip3r_par.ah42)*(RyR_sparse(1+numofryrs:end,2)==2)+ip3r_par.ah42;

                cm=(120*(Yn1(RyR_sparse(ip3rVector),2)/200)-Yn1(RyR_sparse(ip3rVector),1)/10).*(RyR_sparse(1+numofryrs:end,2)==2)+Yn1(RyR_sparse(ip3rVector),1)/10;
                [m42_new,h42_new,m24_new,h24_new]=update_ip3r_states((cm*scale),m42,h42,m24,h24,lambda_h42,ip3r_par,par);

                tan1Vector = [[zeros(numofryrs,1)...
                    min(3.17*1e2*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                    max(0.250*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)];...
                    [zeros(numofip3rs,1)...
                    ip3r_par.a42+ip3r_par.v42.*m42_new.*h42_new...
                    (ip3r_par.a24+ip3r_par.v24.*(1-m24_new.*h24_new)).*ip3r_par.scale...
                    zeros(numofip3rs,1)]];

                tanIndices = [(RyR_sparse(:,2)*2)-1 RyR_sparse(:,2)*2]  ;   

                tanVectorState = tanVector(([tanIndices-1]*total_rec)+[1:total_rec]');
                tan1VectorState = tan1Vector(([tanIndices-1]*total_rec)+[1:total_rec]');

                RyR_sparse(:,3)=RyR_sparse(:,3)+...
                    par.dt.*...
                    (sum(tanVectorState,2)+sum(tan1VectorState,2))/2;



               % Y0=[Cai, Cajsr, Cansr, CaM, CaCaM, ATP, CaATP, F4, F4Ca, TnC];
            %      1     2      3     4     5     6     7     8    9    10  
                Cai=Yn1(:,1);
                Cajsr=Yn1(:,2);
                Cansr=Yn1(:,3);
                CaM=Yn1(:,4);
                CaCaM=Yn1(:,5);
                ATP=Yn1(:,6);
                CaATP=Yn1(:,7);
                F4=Yn1(:,8);
                CaF4=Yn1(:,9);
                TnC=Yn1(:,10);
                h24=h24_new;
                m24=m24_new;
                h42=h42_new;
                m42=m42_new;

                % Determine which state changing receptor should enter
                r1=RyR_sparse(indx_dt_min,1);
                r2=RyR_sparse(indx_dt_min,2);
                Tan1=[tan1Vector(indx_dt_min,1:2);...
                        tan1Vector(indx_dt_min,3:4)];

                Prev_state=Tan1(r2,:);

                new_state=find((cumsum(Prev_state)/sum(Prev_state))>=rand,1); 
    %             new_state=find((cumsum(Prev_state)/sum(Prev_state))>=0.1,1); 
% %                 RyR_sparse(indx_dt_min,2)=new_state;
% %                 RyR_sparse(indx_dt_min,3)=0;
% %                 RyR_sparse(indx_dt_min,4)=rand;
% %                 
%                 if (RyR_open_sums(r1)+Ch_mx(new_state)<0)&& (indx_dt_min<=numofryrs)
%                    disp('problem!') 
%                    
%                 end
                if indx_dt_min<=numofryrs
                    RyR_open_sums(r1)=RyR_open_sums(r1)+Ch_mx(new_state);
                else
                    IP3R_open_sums(r1)=IP3R_open_sums(r1)+Ch_mx(new_state);
                end
                
                                RyR_sparse(indx_dt_min,2)=new_state;
                RyR_sparse(indx_dt_min,3)=0;
                RyR_sparse(indx_dt_min,4)=rand;
                

                end
        dumt=dumt+par.dt;
        end

        % RyR &SR related regions
        RyR_open_save(:,iter)=RyR_open_sums;
        IP3R_open_save(:,iter)=IP3R_open_sums;
        Cai_save(:,iter)=Cai;
        Cajsr_save(:,iter)=Cajsr;
        Cansr_save(:,iter)=Cansr;
        Caf4_save(:,iter)=CaF4;
        ATP_save(:,iter)=ATP;
        CaM_save(:,iter)=CaM;
        TpC_save(:,iter)=TnC;

        time(iter)=simt+dumt;
        simt=time(iter);


    end
    toc            
    parsave(['coupl_rep_uncoupl_mult_clusters_Cannell_Cao_ryrs_',dec2str(numofryrs),'_ip3rs_',dec2str(numofip3rs),'_Tsim_',dec2str(simLength),'smaller_set_test_v',dec2str(fop),'_dist_',dec2str(par.dist),'um_tr13.mat'],time,Cai_save,0.,0.,Caf4_save,RyR_open_save,IP3R_open_save,par,ip3r_par)


end
disp('finally!')
end