ryr_nr=270;
% ip3r_nr=20;
% Number of runs?
% frst=201;
% set_nr=300;
frst=1;
set_nr=100;


ct=2000;
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
% for ip3r_nr=[20,30,40,60]
for ip3r_nr=[20]
for aa=frst:set_nr
    
%     file_name=['coupl_rep_uncoupl_mult_clusters_Cannell_Cao_ryrs_',num2str(ryr_nr),'_ip3rs_',num2str(ip3r_nr),'_Tsim_2000smaller_set_test_v',num2str(aa),'_dist_1-8um_tr13.mat'];
%         file_name=['uncoupl_clusters_Cannell_Cao_ryrs_',num2str(ryr_nr),'_ip3rs_',num2str(ip3r_nr),'_Tsim_2000smaller_set_test_v',num2str(aa),'_dist_2um_tr13_no_st2.mat'];
if ip3r_nr==0
    file_name=['one_d_Cannell_Cao_ryrs_90_ip3rs_0_Tsim_3000smaller_set_test_v',num2str(aa+108),'_tr13.mat'];
else
    file_name=['one_d_Cannell_Cao_ryrs_90_ip3rs_',num2str(ip3r_nr),'_Tsim_3000smaller_set_test_v',num2str(aa),'_tr13.mat'];
end
%     if isfile(file_name)
         
        load(file_name);    
        
        parameters(1,aa-frst+1)=par.gryr;
        parameters(2,aa-frst+1)=par.sc;
        parameters(3,aa-frst+1)=par.scale_sercas;
        parameters(4,aa-frst+1)=par.casr_0;

        linescan_length=size(Caf4_save,1);
        
        x_v=((1:linescan_length)*par.dx-4);
        
        dyad_2=linescan_length/2;
        
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
save(['gryr002_preliminary_tr13_ip3r',num2str(ip3r_nr),'_0_st_one_dyad.mat']);
end
disp('done')
