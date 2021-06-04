ryr_nr=270;
% ip3r_nr=20;
% Number of runs?
% frst=201;
% set_nr=300;
frst=1;
set_nr=100;


ct=3000;
parameters=nan(4,set_nr-frst+1);
% almost_there/

% dyad_amp=cell(3,set_nr);
% dyad_dur=cell(3,set_nr);
% dyad_width=cell(3,set_nr);

coup1_ev_nr=nan(set_nr-frst+1,1);
coup2_ev_nr=nan(set_nr-frst+1,1);
coup1_amp_indx=cell(set_nr-frst+1,1);
coup2_amp_indx=cell(set_nr-frst+1,1);
coup1_max_amp=cell(set_nr-frst+1,1);
coup2_max_amp=cell(set_nr-frst+1,1);

uncp_num_of_ev=nan(set_nr-frst+1,1);
uncp_max_amp=cell(set_nr-frst+1,1);
uncp_amp_index=cell(set_nr-frst+1,1);
uncp_duration=cell(set_nr-frst+1,1);
for ip3r_nr=[20,30,40,60]
for aa=frst:set_nr
    
    file_name=['coupl_rep_uncoupl_mult_clusters_Cannell_Cao_ryrs_',num2str(ryr_nr),'_ip3rs_',num2str(ip3r_nr),'_Tsim_',num2str(ct),'smaller_set_test_v',num2str(aa),'_dist_1-8um_tr13.mat'];
%         file_name=['uncoupl_clusters_Cannell_Cao_ryrs_',num2str(ryr_nr),'_ip3rs_',num2str(ip3r_nr),'_Tsim_2000smaller_set_test_v',num2str(aa),'_dist_2um_tr13_no_st2.mat'];
    
%     if isfile(file_name)
         
        load(file_name);    
        
        parameters(1,aa-frst+1)=par.gryr;
        parameters(2,aa-frst+1)=par.sc;
        parameters(3,aa-frst+1)=par.scale_sercas;
        parameters(4,aa-frst+1)=par.casr_0;

        linescan_length=ceil(8/par.dx);
        
        x_v=((1:linescan_length)*par.dx-4);
        
        dyad_1=(linescan_length/2-2/par.dx);%3
        dyad_2=linescan_length/2;
        dyad_3=linescan_length/2+2/par.dx;
        
        %{
        figure;

        subplot(4,2,1)
        imagesc(Caf4_save/Caf4_save(1,1)-1)
        
        subplot(4,2,3)
        plot(mean(Caf4_save(dyad_1-2:dyad_1+2,:))/Caf4_save(1,1)-1)
        
        subplot(4,2,4)
        
        [d1_amp,d1_m]=max(mean(Caf4_save(dyad_1-2:dyad_1+2,:))/Caf4_save(1,1)-1);
        plot(x_v,Caf4_save(:,d1_m)/Caf4_save(1,1)-1)
        hold on
        plot([-4,4],[d1_amp/2,d1_amp/2],'r')
        hold off
        xlim([-4,4])
        
        subplot(4,2,5)
        plot(mean(Caf4_save(dyad_2-2:dyad_2+2,:))/Caf4_save(1,1)-1)
        
        subplot(4,2,6)
        
        [d2_amp,d2_m]=max(mean(Caf4_save(dyad_2-2:dyad_2+2,:))/Caf4_save(1,1)-1);
        

        plot(x_v,Caf4_save(:,d2_m)/Caf4_save(1,1)-1)
        hold on
        plot([-4,4],[d2_amp/2,d2_amp/2],'r')
        hold off
        xlim([-4,4])
        
        subplot(4,2,7)
        plot(mean(Caf4_save(dyad_3-2:dyad_3+2,:))/Caf4_save(1,1)-1)
        
        subplot(4,2,8)
        
        [d3_amp,d3_m]=max(mean(Caf4_save(dyad_3-2:dyad_3+2,:))/Caf4_save(1,1)-1);
        plot(x_v,Caf4_save(:,d3_m)/Caf4_save(1,1)-1)
        hold on
        plot([-4,4],[d3_amp/2,d3_amp/2],'r')
        hold off
        xlim([-4,4])
        
        %}
        
        coup_dyad1=mean(Caf4_save(dyad_1-2:dyad_1+2,:))/Caf4_save(1,1)-1;
        bool_c1=(coup_dyad1>=1/2);
        beg_c1=find(diff(bool_c1)>0);
        end_c1=find(diff(bool_c1)<0);
        if numel(end_c1)<numel(beg_c1)
            end_c1=[end_c1,numel(coup_dyad1)];
        end
        coup1_ev_nr(aa-frst+1,1)=numel(beg_c1);
        
        if coup1_ev_nr(aa-frst+1,1)>0
           
            coup1_amp_indx{aa-frst+1,1}=nan(coup1_ev_nr(aa-frst+1,1),1);
            coup1_max_amp{aa-frst+1,1}=nan(coup1_ev_nr(aa-frst+1,1),1);
            
            for c1=1:coup1_ev_nr(aa-frst+1,1)
               [c1_amp,c1_indx]=max(coup_dyad1(beg_c1(c1):end_c1(c1)));
               coup1_amp_indx{aa-frst+1,1}(c1)=c1_indx+beg_c1(c1);
               coup1_max_amp{aa-frst+1,1}(c1)=c1_amp;
            end
            
        end
        
        coup_dyad2=mean(Caf4_save(dyad_3-2:dyad_3+2,:))/Caf4_save(1,1)-1;
        bool_c2=(coup_dyad2>=1/2);
        beg_c2=find(diff(bool_c2)>0);
        end_c2=find(diff(bool_c2)<0);
        
        if numel(end_c2)<numel(beg_c2)
            end_c2=[end_c2,numel(coup_dyad2)];
        end
        
       coup2_ev_nr(aa-frst+1,1)=numel(beg_c2);
        
        if coup2_ev_nr(aa-frst+1,1)>0
           
            coup2_amp_indx{aa-frst+1,1}=nan(coup2_ev_nr(aa-frst+1,1),1);
            coup2_max_amp{aa-frst+1,1}=nan(coup2_ev_nr(aa-frst+1,1),1);
            
            for c2=1:coup2_ev_nr(aa-frst+1,1)
               [c2_amp,c2_indx]=max(coup_dyad2(beg_c2(c2):end_c2(c2)));
               coup2_amp_indx{aa-frst+1,1}(c2)=c2_indx+beg_c2(c2);
               coup2_max_amp{aa-frst+1,1}(c2)=c2_amp;
            end
            
        end
        % Look for events - times when ca amplitude exceeds 1/2
        uncoupled_dyad=mean(Caf4_save(dyad_2-2:dyad_2+2,:))/Caf4_save(1,1)-1;
        bool_array=(uncoupled_dyad>=1/2);
        dummy2=diff(bool_array);
        beg_ev=find(dummy2>0);
        end_ev=find(dummy2<0);
        uncp_num_of_ev(aa-frst+1)=numel(beg_ev);
        full_ev=min(numel(beg_ev),numel(end_ev));
        % If there are any events, look into them
        if full_ev>0
            uncp_max_amp{aa-frst+1}=nan(full_ev,1);
            uncp_amp_index{aa-frst+1}=nan(full_ev,1);
            uncp_duration{aa-frst+1}=nan(full_ev,1);
            
            indx1=0;
            for bb=1:full_ev
                % Find the amplitude of the event
                [unc_amp,unc_indx]=max(uncoupled_dyad(beg_ev(bb):end_ev(bb)));
                uncp_max_amp{aa-frst+1}(bb)=unc_amp;
                uncp_amp_index{aa-frst+1}(bb)=beg_ev(bb)+unc_indx;
                % Added min(1/2, unc_amp) bc large events were prematurely
                % ending
%                 uncp_bool_array=(uncoupled_dyad(indx1+1:end)>=max(min(unc_amp/2,1/2),min(uncoupled_dyad(indx1+1:end_ev(bb)))+1e-4));
                uncp_bool_array=(uncoupled_dyad(indx1+1:end)>=unc_amp/2);
                uncp_dm=diff(uncp_bool_array);
                in1=find(uncp_dm>0,1,'first');
                in2=find(uncp_dm(in1:end)<0,1,'first');
                if isempty(in1)|isempty(in2)
                    continue
                end
                uncp_duration{aa-frst+1}(bb)=time(in2+in1+indx1)-time(in1+indx1);
                indx1=indx1+in2+in1;
            end
        
        end
        
%     end
end
save(['gryr002_preliminary_tr13_ip3r',num2str(ip3r_nr),'_3_st_all_dyads.mat']);
end
disp('done')
