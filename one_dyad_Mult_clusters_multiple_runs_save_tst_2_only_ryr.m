% Agne Tilunaite 2019
% Edited for clarity by Hilary Hunt 2020
numSimulations=100;
simLength=3000; %ms

IP3Conc=0.15; %uM Note - should be unnecessary
casr0=800*ones(numSimulations,1);
num_ip3r=0;
parfor(fop=1:numel(casr0),20)
                
    par=sobie_etal_2002_combined_parameters_3(0.1,casr0(fop));
    par.casr_0=casr0(fop);
%     par.sc=examine_sc(fop);
    par.sc=10;
%     par.scale_sercas=ex_serca(fop);
    par.scale_sercas=1;
    par.gillc=0;
%     par.gillc=1.5e-16;
%     par.gillc=gilcc(fop);
%     par.gryr=0.015;
%      par.gryr=ex_ryr(fop);
    par.gryr=0.02;
    par.dt0=1e-4; %ms
    obsT=simLength; %ms
    % linesca6length=150;
    par.dx=0.2/5; %um
  
    par.num_of_varied_ip3r=num_ip3r;
    
    par.ryr_op_times=800;
    par.keep_ryr_open=5;
   
    c0=0.1;
    % Added blank
    ip3r_par=Cao_et_al_2014(1,c0,IP3Conc);
    ip3r_par.kipr=0.02;
%     par.kip3r=ip3r_par.kipr;

    linescan_length=ceil(2/par.dx);

    FOV=[1,linescan_length];

    RyR_place=ones(FOV(1),FOV(2)); % In other words RyR_place=RyR_matrix>0

    % Storage of each receptor (RyR and then IP3R) behaviour and its transition times:
    % Flatten information about all receptors:
    % 1st column - receptor position in the flattened receptor matrix
    % 2nd column - the state of the receptor
    % 3rd column - waiting time till receptor changes state
    % 4th column - random number to determine receptor transition times

    numofryrs=90;
    numofip3rs=2*par.num_of_varied_ip3r;

    total_rec=numofryrs+numofip3rs;


    RyR_release_place=ceil(linescan_length/2); % for uncoupled one
    RyR_sparse=zeros(numofryrs+numofip3rs,4);
    % RyR_pos=RyR_release_place*ones(numofryrs,1);

    RyR_open_sums=zeros(linescan_length,1);
    IP3R_open_sums=0;

    initiationT=obsT;
    initiationiter=ceil(initiationT/par.dt0);
    titer=ceil(obsT/par.dt0)/100;

    % indx_fl=0;
    % for aa=1:numel(RyR_pos)
    %    RyR_sparse(indx_fl+1:indx_fl+1,1)=RyR_pos(aa);
    %    indx_fl=indx_fl+1;
    % end


    
    div=numofryrs/2;

    RyR_sparse(1:div,1)=RyR_release_place-1;       %1st coupled cluster
    RyR_sparse(div+1:2*div,1)=RyR_release_place+1;


    par.positions=zeros(linescan_length,1);
    par.positions(RyR_release_place-2:RyR_release_place+2)=1;

    par.sercas=~par.positions;
     
    %State 1 closed; state 2 open
    Ch_mx=[-1;1];

    RyR_sparse(:,2)=1; % state of each receptor
    RyR_sparse(1:2,2)=1; % state of each receptor
    RyR_sparse(div+1:div+2,2)=1; % state of each receptor

% Comment if no trigger is applied:
%     RyR_sparse(1:30,2)=2; % state of each receptor
%     RyR_sparse(div+1:div+30,2)=2; % state of each receptor

    
    
    RyR_open_sums(RyR_sparse(1,1))=sum(RyR_sparse(1:div,2)==2);
    RyR_open_sums(RyR_sparse(div+1,1))=sum(RyR_sparse(div+1:2*div,2)==2);
    


    RyR_sparse(:,4)=rand(total_rec,1);

    % Variables I will need to track

    % Y0=[Cai, Cajsr, Cansr, CaM, CaCaM, ATP, CaATP, F4, F4Ca, TnC];
    %      1     2      3     4     5     6     7     8    9    10   

    % RyR &SR related regions

   

    RyR_open_save=zeros(linescan_length,titer);
    IP3R_open_save=0;
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

    tic
    par.dt=par.dt0;
    simt=0;
    rep=1;
    old_indx=0;
    for iter=1:titer
        dumt=0;
        while dumt<=(par.dt0*100)
           
            if (simt+dumt)>(rep*par.ryr_op_times) && (simt+dumt)<(rep*par.ryr_op_times+par.keep_ryr_open)
%                 RyR_sparse(1:30,2)=2; % state of each receptor
%                 RyR_sparse(div+1:div+30,2)=2; % state of each receptor
      
                %RyR_sparse(1:2,2)=2; % state of each receptor
                %RyR_sparse(div+1:div+2,2)=2; % state of each receptor

                RyR_open_sums(RyR_sparse(1,1))=sum(RyR_sparse(1:div,2)==2);
                RyR_open_sums(RyR_sparse(div+1,1))=sum(RyR_sparse(div+1:2*div,2)==2);
                
                rep=ceil((simt+dumt)/par.ryr_op_times);

            end
            
           Yn=[Cai, Cajsr, Cansr, CaM,CaCaM, ATP, CaATP, F4, CaF4, TnC];
           Yn1=RK4_Sobie_Cao_2_diffusion(Yn,par,ip3r_par,RyR_open_sums,0,simt);             % values at the time step n+1


            % Dumm0 variable to find out real minimal time step between the state
            % change of each receptor.
% %             dt1=par.dt*ones(total_rec,1);
            dt1=par.dt0*ones(total_rec,1);



            %vectorised code again
            ryrVector = (1:numofryrs);
            tanVector = [[zeros(numofryrs,1)...
                min(3.17*1e2*((Yn(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                max(0.250*((Yn(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)]];


            tan1Vector = [[zeros(numofryrs,1)...
                min(3.17*1e2*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                max(0.250*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)]];

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
                Yn1=RK4_Sobie_Cao_2_diffusion(Yn,par,ip3r_par,RyR_open_sums,0,simt);          % values at the time step n+1      
                if indx_dt_min==0
                    par.dt=par.dt0;
                    
                else
                   
%                     disp([dt_min*1000,indx_dt_min,par.dt*1000])
                    
    %             %replace the loop above with vectorised code
                ryrVector = (1:numofryrs);
                
                tanVector = [[zeros(numofryrs,1)...
                    min(3.17*1e2*((Yn(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                    max(0.250*((Yn(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)]];

                tan1Vector = [[zeros(numofryrs,1)...
                    min(3.17*1e2*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^2.8),0.7) ...
                    max(0.250*((Yn1(RyR_sparse(ryrVector),1)/par.sc).^-0.5),0.9) zeros(numofryrs,1)]];

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


                % Determine which state changing receptor should enter
                r1=RyR_sparse(indx_dt_min,1);
                r2=RyR_sparse(indx_dt_min,2);
                Tan1=[tan1Vector(indx_dt_min,1:2);...
                        tan1Vector(indx_dt_min,3:4)];

                Prev_state=Tan1(r2,:);

                new_state=find((cumsum(Prev_state)/sum(Prev_state))>=rand,1); 
                
                RyR_open_sums(r1)=RyR_open_sums(r1)+Ch_mx(new_state);
                RyR_sparse(indx_dt_min,2)=new_state;
                RyR_sparse(indx_dt_min,3)=0;
                RyR_sparse(indx_dt_min,4)=rand;
                

                end
        dumt=dumt+par.dt;
        end

        % RyR &SR related regions
        RyR_open_save(:,iter)=RyR_open_sums;
        IP3R_open_save=0;
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
%     parsave(['one_d_Cannell_Cao_ryrs_',num2str(numofryrs),'_ip3rs_',num2str(numofip3rs),'_Tsim_',num2str(obsT),'smaller_set_test_v',num2str(fop),'_tr13_no_trigger_at_all.mat'],time,Cai_save,Cajsr_save,Cansr_save,Caf4_save,RyR_open_save,IP3R_open_save,par)
% Added blank
%     parsave(['one_d_Cannell_Cao_ryrs_',num2str(numofryrs),'_ip3rs_',num2str(numofip3rs),'_Tsim_',num2str(obsT),'smaller_set_test_v',num2str(fop),'_tr13_try_initiate_2ryrs.mat'],time,Cai_save,Cajsr_save,Cansr_save,Caf4_save,RyR_open_save,IP3R_open_save,par,[])
parsave(['one_d_Cannell_Cao_ryrs_',num2str(numofryrs),'_ip3rs_0_Tsim_',num2str(obsT),'smaller_set_test_v',num2str(fop),'_tr13.mat'],time,Cai_save,Cajsr_save,Cansr_save,Caf4_save,RyR_open_save,IP3R_open_save,par,ip3r_par)

end
disp('finally!')
