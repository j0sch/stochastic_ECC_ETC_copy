
% Point to files
% addpath('/Volumes/UniWD/Agne/ip3sparks/')
% addpath('../sims/')
addpath('./helping_functions/')
addpath('/Volumes/UniWD/Agne/ip3sparks/sims_gryr/')
ip3r=[20,30,40,60];

cn=numel(ip3r);


near_distd=cell(numel(ip3r),4);

st_t={'_3','_0','_0'};
runs_with_sparks=zeros(3,numel(ip3r));
for ppp=1:3
for cc=1:numel(ip3r)
% Load condensed simulation files
if ppp~=3   
    load(['gryr002_preliminary_tr13_ip3r',num2str(ip3r(cc)),st_t{ppp},'_st_all_dyads.mat']);
else
    ip3r=[0,10,20,40];
    load(['gryr002_preliminary_tr13_ip3r',num2str(ip3r(cc)),st_t{ppp},'_st_one_dyad.mat']);
    ip3r=[20,30,40,60];
end
    
run_nr=set_nr-frst+1;

near_dist_d=cell(run_nr,1);

for aa=1:run_nr
    if numel(uncp_amp_index{aa})>0
        near_dist_d{aa,1}=nan(numel(uncp_amp_index{aa}),1);

            for bb=1: numel(uncp_amp_index{aa})
               dm1=find((coup1_amp_indx{aa}-uncp_amp_index{aa}(bb))<=0,1,'last');
               dm2=find((coup2_amp_indx{aa}-uncp_amp_index{aa}(bb))<=0,1,'last'); 
               if isempty(dm1)&& isempty(dm2)
                 near_dist_d{aa,1}(bb,1)=time(uncp_amp_index{aa}(bb));
               elseif isempty(dm2)&&(isempty(dm1)==0)
                 near_dist_d{aa,1}(bb,1)=(-time(coup1_amp_indx{aa}(dm1,1))+time(uncp_amp_index{aa}(bb)));
               elseif isempty(dm1) && isempty(dm2)==0
                 near_dist_d{aa,1}(bb,1)=(-time(coup2_amp_indx{aa}(dm2,1))+time(uncp_amp_index{aa}(bb)));
               else
                near_dist_d{aa,1}(bb,1)=min(-time(coup1_amp_indx{aa}(dm1,1))+time(uncp_amp_index{aa}(bb)),...
                   -time(coup2_amp_indx{aa}(dm2,1))+time(uncp_amp_index{aa}(bb)));
               end
            end
    end 
end
runs_with_sparks(ppp,cc)=100-sum(cellfun(@(x)isempty(x),near_dist_d));
near_distd{cc,ppp}=cell2mat(near_dist_d);
end
end

if 0
run_nr=500;
rand_sampling=nan(run_nr,3);
for aa=1:run_nr
    rand_sampling(aa,:)=randperm(100,3);
end

for cc=1:numel(ip3r)
    
    load(['gryr002_preliminary_tr13_ip3r',num2str(ip3r(cc)-20+2),'_no_st_all_dyads_dummy_dyads_1.mat'],'uncp_amp_index','time');
    
    near_distd{cc,4}=time(cell2mat(uncp_amp_index));
    
    
    near_dist_d=cell(run_nr,1);


    for aa=1:run_nr
        if numel(uncp_amp_index{rand_sampling(aa,2)})>0
            near_dist_d{aa,1}=nan(numel(uncp_amp_index{rand_sampling(aa,2)}),1);

                for bb=1: numel(uncp_amp_index{rand_sampling(aa,2)})
                   dm1=find((uncp_amp_index{rand_sampling(aa,1)}-uncp_amp_index{rand_sampling(aa,2)}(bb))<=0,1,'last');
                   dm2=find((uncp_amp_index{rand_sampling(aa,3)}-uncp_amp_index{rand_sampling(aa,2)}(bb))<=0,1,'last'); 
                   if isempty(dm1)&& isempty(dm2)
                     near_dist_d{aa,1}(bb,1)=time(uncp_amp_index{rand_sampling(aa,2)}(bb));
                   elseif isempty(dm2)&&(isempty(dm1)==0)
                     near_dist_d{aa,1}(bb,1)=(-time(uncp_amp_index{rand_sampling(aa,1)}(dm1,1))+time(uncp_amp_index{rand_sampling(aa,2)}(bb)));
                   elseif isempty(dm1) && isempty(dm2)==0
                     near_dist_d{aa,1}(bb,1)=(-time(uncp_amp_index{rand_sampling(aa,3)}(dm2,1))+time(uncp_amp_index{rand_sampling(aa,2)}(bb)));
                   else
                    near_dist_d{aa,1}(bb,1)=min(-time(uncp_amp_index{rand_sampling(aa,1)}(dm1,1))+time(uncp_amp_index{rand_sampling(aa,2)}(bb)),...
                       -time(uncp_amp_index{rand_sampling(aa,3)}(dm2,1))+time(uncp_amp_index{rand_sampling(aa,2)}(bb)));
                   end
                end
        end 
    end

near_distd{cc,3}=cell2mat(near_dist_d);
    
end
end
%%

% dmm=reshape(near_distd,[1 4]);
% 
% plotSpread(near_distd,'xNames',{'ip3r=0','ip3r=5','10','ip3r=20'})
% hold on
% violin(dmm,'bw',100)
% hold off
% ylim([0,2000])


%%

figure;
for ccc=1:cn
barnr=20;
intw=2000/barnr;
bars=zeros(barnr,4);

for aaa=1:barnr
   for bbb=1:4
      bars(aaa,bbb)=sum((near_distd{ccc,bbb}>=(aaa-1)*intw).*(near_distd{ccc,bbb}<aaa*intw))/numel(near_distd{ccc,bbb}); 
   end
end
 subplot(1,cn,ccc)
 bar(bars)
 title(['ip3r=',num2str((ip3r(ccc)-20)/2)])
end

%%
figure;
for ccc=1:cn
% barnr=20;
% intw=2000/barnr;
barnr=5;
intw=50;
bars=zeros(barnr,4);

for aaa=1:barnr
   for bbb=1:4
      bars(aaa,bbb)=(sum((near_distd{ccc,bbb}<=(aaa*intw)))/numel(near_distd{ccc,bbb}))*100; 
   end
end
 subplot(1,cn,ccc)
 hbar1=bar(bars);
 ylabel('% of ISI')
 set(gca, 'XTickLabel',{['t<=',num2str(intw)],['t<=',num2str(2*intw)],...
     ['t<=',num2str(3*intw)],['t<=',num2str(4*intw)],['t<=',num2str(5*intw)]})
 ylim([0,100])
 xtickangle(30)
 title(['ip3r=',num2str((ip3r(ccc)-20)/2)])
end
fcolor=cell2mat(get(hbar1,'FaceColor'));

% fcolor=fcolor2;
%%
figure;

cn=numel(ip3r);
plot_vio=near_distd(:,cc);
for cc=1:cn
%     subplot(1,cn,cc)
subplot(1,cn,cc)
plotSpread(near_distd(cc,1),'xNames',{'3tr 3 dyads','0 tr 3 dyads','Shuffled control','0 tr 1 dyad'})
hold on
violin(near_distd(cc,1),'facecolor',fcolor,'bw',100,'mc',[],'medc',[])
% plot([0,4],[intw*barnr, intw*barnr],'b--','Linewidth',3)
hold off
ylim([0,2000])
ylabel('ISI (ms)')
xtickangle(30)
% ylabel('time between sparks in different regions (ms)')
title(['IP_3R',{' '},num2str((ip3r(cc)-20)/2)])
end

%%
figure;

cn=numel(ip3r);
% Edited - near_distd(c,:) -> near_distd(cc,1)
for cc=1:cn
%     subplot(1,cn,cc)
subplot(1,cn,cc)
plotSpread(near_distd(cc,1),'xNames',{''})
hold on
violin(near_distd(cc,1),'facecolor',fcolor,'bw',100,'mc',[],'medc',[])
% plot([0,4],[intw*barnr, intw*barnr],'b--','Linewidth',3)
hold off
ylim([0,600])
ylabel('ISI (ms)')
xtickangle(30)
% ylabel('time between sparks in different regions (ms)')
title(['ip3r=',num2str((ip3r(cc)-20)/2)])
end
%% Wilcoxon rank sum pairwise test

pv=nan(cn,3);
hv=nan(cn,3);

for cc=1:cn
   [pv(cc,1),hv(cc,1)]=ranksum(near_distd{cc,1},near_distd{cc,2});
   [pv(cc,2),hv(cc,2)]=ranksum(near_distd{cc,1},near_distd{cc,3});
   [pv(cc,3),hv(cc,3)]=ranksum(near_distd{cc,2},near_distd{cc,3});
end


%% Horizontal column laying. Scatter and cumulative probability together

figure;
for ccc=1:cn
% barnr=20;
% intw=2000/barnr;
barnr=5;
intw=50;
bars=zeros(barnr,4);

for aaa=1:barnr
   for bbb=1:4
      bars(aaa,bbb)=(sum((near_distd{ccc,bbb}<=(aaa*intw)))/numel(near_distd{ccc,bbb}))*100; 
   end
end
 subplot(4,2,2*ccc)
 hbar1=bar(bars);
 ylabel('% of ISI')
 set(gca, 'XTickLabel',{['t<=',num2str(intw)],['t<=',num2str(2*intw)],...
     ['t<=',num2str(3*intw)],['t<=',num2str(4*intw)],['t<=',num2str(5*intw)]})
 
 ylim([0,100])
 xtickangle(15)
%  title(['ip3r=',num2str((ip3r(ccc)-20)/2)])
if ccc==1
    title('Percentage of ISI <= t')
end
 ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end
fcolor=cell2mat(get(hbar1,'FaceColor'));
% subplot(2,4,5)
% ylabel('P(spark in t ms)')


cn=numel(ip3r);

for cc=1:cn
%     subplot(1,cn,cc)
subplot(4,2,2*cc-1)
plotSpread(near_distd(cc,:),'xNames',{'triggered','untriggered','control','1 dyad'})
hold on
violin(near_distd(cc,:),'facecolor',fcolor,'bw',100,'mc',[],'medc',[])
% plot([0,4],[intw*barnr, intw*barnr],'b--','Linewidth',3)
hold off
xtickangle(15)
ylim([0,2000])
% ylabel('Inter spark interval (ms)')
ylabel('ISI (ms) ')
% ylabel('time between sparks(ms)')
% ylabel('time between sparks in different regions (ms)')
% title(['ip3r=',num2str((ip3r(cc)-20)/2)])
if cc==1
    title('Time between sparks') 
end
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end
% subplot(2,4,1)
% ylabel('time between sparks(ms)')

%% Bar chart

%% Horizontal column laying. Scatter and cumulative probability together

figure;
for ccc=1:cn
% barnr=20;
% intw=2000/barnr;
barnr=5;
intw=50;
bars=zeros(barnr,4);

for aaa=1:barnr
   for bbb=1:4
      bars(aaa,bbb)=(sum((near_distd{ccc,bbb}<=(aaa*intw)))/numel(near_distd{ccc,bbb}))*100; 
   end
end
 subplot(4,1,ccc)
 hbar1=bar(bars);
 ylabel('% of ISI')
 set(gca, 'XTickLabel',{['t<=',num2str(intw)],['t<=',num2str(2*intw)],...
     ['t<=',num2str(3*intw)],['t<=',num2str(4*intw)],['t<=',num2str(5*intw)]})
 
 ylim([0,100])
 xtickangle(15)
%  title(['ip3r=',num2str((ip3r(ccc)-20)/2)])
if ccc==1
    title('Percentage of ISI <= t')
end
 ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end
fcolor=cell2mat(get(hbar1,'FaceColor'));
% subplot(2,4,5)
% ylabel('P(spark in t ms)')


%% Bar chart for Llew's paper
xlabs=cell(1,4);
for ccc=1:4
xlabs{ccc}=['# IP_3Rs =',' ',num2str((ip3r(ccc)-20)/2)];
end
xlabc=categorical(xlabs);
xlabc=reordercats(xlabc,xlabs);
figure
bar(xlabc,runs_with_sparks',0.8)
ylabel('% of runs with sparks')
set(gca,'FontSize',16)

%% Heat map of sim - use to visualise individual sims
figure('Position',[10,10,700,350])
surf(time/1000,(1:size(Caf4_save,1))*0.04,Caf4_save,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');view(0,90)
xlabel('time (s)')
ylabel('distance (\mum)')
axis('tight')
% colorbar
% ylabel('CaF4 (\mu M)')
set(gca,'FontSize',16)