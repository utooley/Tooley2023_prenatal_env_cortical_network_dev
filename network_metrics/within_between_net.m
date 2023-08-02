%% Do network metrics for birth data
datadir=fullfile('~/data/smyser/smyser1/wunder/eLABe/gordon_pconns_plus_atlas_subcortical/full_mats/')
outdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/birth/'

subjlistdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/subjLists/'
subjlist=readtable(fullfile(subjlistdir,'n319_birth_fc_ids.csv'),'Delimiter',',','ReadVariableNames', 1)

% these pconns are already z-scored
clear modul
clear avgweight
clear num_comms_modul
clear system_segreg
clear mean_within_sys
clear mean_between_sys
clear system_conn
clear part_coef_pos
clear part_coef_neg
clear avgclustco_both
clear avgclustco_all
clear part_coef_avg_all
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);
part_coef_pos=zeros(height(subjlist),1);
part_coef_neg=zeros(height(subjlist),1);
avgclustco_both=zeros(height(subjlist),1);
avgclustco_all=zeros(height(subjlist),333);
part_coef_avg_all=zeros(height(subjlist),333);
gordon_parcels=readtable(fullfile('~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.network.csv'));
%% Load FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.x(n)) %look at this
    file=fullfile(datadir,strcat(sub,'_gordon_parcel_plus_term_N50_eLABe_atlas_subcort.pconn.nii'));
    subfcmat=cifti_read(file).cdata; 
    %eliminate the 7.25s on the diagonal
    for x=1:size(subfcmat)
       subfcmat(x,x)=0;
    end
    %make it 333 X 333
    subfcmat=subfcmat(1:333,1:333);
    avgweight(n,1)=mean(subfcmat(subfcmat~=0));
    
    %segregation proportionally weighted
    [S, W, B] = segregation(subfcmat,gordon_parcels.numeric_network);
    system_segreg_prop(n,1)=S;
    mean_within_sys_prop(n,1)=W;
    mean_between_sys_prop(n,1)=B;

    %Within and between connectivity, code adapted from Micaela Chan 2018
    Ci=gordon_parcels.numeric_network;
    nCi = unique(Ci);
    M=subfcmat;
    
    for i = 1:length(nCi) % loop through communities
        for j = 1:length(nCi)
           Wi = Ci == nCi(i); % find index for within community edges
           Bi = Ci == nCi(j); % find index for between community edges to specific community
           
           Wv_temp = M(Wi,Wi); % extract within community edges
           Bv_temp = M(Wi,Bi); % extract between community edges to specific community
           
           %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
           Bv = [Bv_temp(:)'];
           system_connectivity(i,j)=mean(Bv(Bv~=0));
           %system_between(i,1)=mean(Bv);
           %if i==j
           %else
        end
    end

%transpose these (or something so they can be saved out on a subject basis.
system_connectivity;
%put this all into a matrix for everyone so that we can see the average
%system connectivity
system_connectivity_all(:,:,n)=system_connectivity;
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

%modularity using Gen Louvain classic modularity maximization, negative
%asymmetric weighting
for c = 1:100 %run it 100x for each subject, get an average number of communities detected.
    [M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
    modul_temp(c)=Q;
    num_comms_temp(c)=length(unique(M));
%also save the number of communities detected.
end
modul(n,1)=mean(modul_temp(:));
num_comms_modul(n,1)=mean(num_comms_temp(:));
%average participation coefficient!
[Ppos, Pneg]=participation_coef_sign(subfcmat, gordon_parcels.numeric_network);
part_coef_pos(n,1)=mean(Ppos);
part_coef_neg(n,1)=mean(Pneg);
%write out all nodes participation coefficient
part_coef_avg_all(n,:)=((Ppos+Pneg)/2);

%average clustering coefficient
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));
%write out all nodes clustering coefficient
avgclustco_all(n,:)=clustering_coef_wu_sign(subfcmat,3);
end

save('~/birth_network_metrics_333_n319.mat')

header={'ID', 'avgweight', 'modul_avg', 'avgclustco_both','num_comms_modul_avg','part_coef_pos','part_coef_neg', 'system_segreg', 'mean_within_sys', 'mean_between_sys', compose('net%d', 1:78)}

outfile=table(char(subjlist.x), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg_prop, mean_within_sys_prop, mean_between_sys_prop, system_conn)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, 'n319_nreg333_birth_within_between_gordon_withmodulpartcoef.mat'), 'outfile2')
writetable(outfile2,fullfile(outdir,'n319_nreg333_birth_within_between_gordon_withmodulpartcoef.csv'))

%also save the mean system connectivity matrix just in case
mean_system_conn_mat=mean(system_connectivity_all,3)
save(fullfile(outdir, strcat('n319_nreg333_birth_mean_system_conn.mat')), 'mean_system_conn_mat')

%save the nodal participation coefficient and clustering coefficient
outfile=table(char(subjlist.x), part_coef_avg_all)
writetable(outfile,strcat(outdir,'/n319_nreg333_birth_part_coef_avg_nodewise.csv'))

outfile=table(char(subjlist.x), avgclustco_all)
writetable(outfile,strcat(outdir,'/n319_nreg333_birth_clust_coef_avg_nodewise.csv'))

%% Do network metrics for year 2 data

datadir=fullfile('~/smyser3/wunder/eLABe_II_fc/gordon_ptseries_pconns/subcort_pconns/')
outdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/y2/'

subjlistdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/subjLists/'
subjlist=readtable(fullfile(subjlistdir,'n114_y2_fc_ids.csv'),'Delimiter',',','ReadVariableNames', 1)

% these pconns are already z-scored
clear modul
clear avgweight
clear num_comms_modul
clear system_segreg
clear mean_within_sys
clear mean_between_sys
clear system_conn
clear part_coef_pos
clear part_coef_neg
clear avgclustco_both
clear avgclustco_all
clear part_coef_avg_all
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);
part_coef_pos=zeros(height(subjlist),1);
part_coef_neg=zeros(height(subjlist),1);
avgclustco_both=zeros(height(subjlist),1);
avgclustco_all=zeros(height(subjlist),333);
part_coef_avg_all=zeros(height(subjlist),333);
gordon_parcels=readtable(fullfile('~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.network.csv'));
%% Load FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.long_id(n)) %look at this
    file=fullfile(datadir,strcat(sub,'_gordon_parcels_plus_subcort.pconn.nii'));
    subfcmat=cifti_read(file).cdata;
    %eliminate the 7.25s on the diagonal
    for x=1:size(subfcmat)
       subfcmat(x,x)=0;
    end
    %make it 333 X 333
    subfcmat=subfcmat(1:333,1:333);
    avgweight(n,1)=mean(subfcmat(subfcmat~=0));
    
    %segregation proportionally weighted
    [S, W, B] = segregation(subfcmat,gordon_parcels.numeric_network);
    system_segreg_prop(n,1)=S;
    mean_within_sys_prop(n,1)=W;
    mean_between_sys_prop(n,1)=B;

    %Within and between connectivity, code adapted from Micaela Chan 2018
    Ci=gordon_parcels.numeric_network;
    nCi = unique(Ci);
    M=subfcmat;
    
    for i = 1:length(nCi) % loop through communities
        for j = 1:length(nCi)
           Wi = Ci == nCi(i); % find index for within community edges
           Bi = Ci == nCi(j); % find index for between community edges to specific community
           
           Wv_temp = M(Wi,Wi); % extract within community edges
           Bv_temp = M(Wi,Bi); % extract between community edges to specific community
           
           %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
           Bv = [Bv_temp(:)'];
           system_connectivity(i,j)=mean(Bv(Bv~=0));
           %system_between(i,1)=mean(Bv);
           %if i==j
           %else
        end
    end

%transpose these (or something so they can be saved out on a subject basis.
system_connectivity;
%put this all into a matrix for everyone so that we can see the average
%system connectivity
system_connectivity_all(:,:,n)=system_connectivity;
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

%modularity using Gen Louvain classic modularity maximization, negative
%asymmetric weighting
for c = 1:100 %run it 100x for each subject, get an average number of communities detected.
    [M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
    modul_temp(c)=Q;
    num_comms_temp(c)=length(unique(M));
%also save the number of communities detected.
end
modul(n,1)=mean(modul_temp(:));
num_comms_modul(n,1)=mean(num_comms_temp(:));
%average participation coefficient!
[Ppos, Pneg]=participation_coef_sign(subfcmat, gordon_parcels.numeric_network);
part_coef_pos(n,1)=mean(Ppos);
part_coef_neg(n,1)=mean(Pneg);
%write out all nodes participation coefficient
part_coef_avg_all(n,:)=((Ppos+Pneg)/2);

%average clustering coefficient
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));
%write out all nodes clustering coefficient
avgclustco_all(n,:)=clustering_coef_wu_sign(subfcmat,3);
end

save('~/Box/projects/in_progress/n114_y2_network_metrics.mat')

header={'ID', 'avgweight', 'modul_avg', 'avg_clustco_both','num_comms_modul_avg','part_coef_pos','part_coef_neg', 'system_segreg', 'mean_within_sys', 'mean_between_sys', compose('net%d', 1:78)}

outfile=table(char(subjlist.modid), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg_prop, mean_within_sys_prop, mean_between_sys_prop, system_conn)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, 'n114_y2_within_between_gordon_withmodulpartcoef.mat'), 'outfile2')
writetable(outfile2,fullfile(outdir,'n114_y2_within_between_gordon_withmodulpartcoef.csv'))

%also save the mean system connectivity matrix
mean_system_conn_mat=mean(system_connectivity_all,3)
save(fullfile(outdir, strcat('n114_y2_mean_system_conn.mat')), 'mean_system_conn_mat')

%save the nodal participation coefficient and clustering coefficient
outfile=table(char(subjlist.modid), part_coef_avg_all)
writetable(outfile,strcat(outdir,'/n114_y2_part_coef_avg_nodewise.csv'))

outfile=table(char(subjlist.modid), avgclustco_all)
writetable(outfile,strcat(outdir,'/n114_y2_clust_coef_avg_nodewise.csv'))

%% Do network metrics for year 3 data
datadir=fullfile('~/10.20.145.4/SMYSER03/smyser3/wunder/eLABe_III_fc/gordon_ptseries_pconns/subcort_pconns/')
outdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/y3/'

subjlistdir='~/Box/projects/in_progress/within_between_network_longitudinal/data/subjLists/'
subjlist=readtable(fullfile(subjlistdir,'n88_y3_fc_ids.csv'),'Delimiter',',','ReadVariableNames', 1)

% these pconns are already z-scored
clear modul
clear avgweight
clear num_comms_modul
clear system_segreg
clear mean_within_sys
clear mean_between_sys
clear system_conn
clear part_coef_pos
clear part_coef_neg
clear avgclustco_both
clear avgclustco_all
clear part_coef_avg_all
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);
part_coef_pos=zeros(height(subjlist),1);
part_coef_neg=zeros(height(subjlist),1);
avgclustco_both=zeros(height(subjlist),1);
avgclustco_all=zeros(height(subjlist),333);
part_coef_avg_all=zeros(height(subjlist),333);
gordon_parcels=readtable(fullfile('~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.network.csv'));
%% Load FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.modid(n)) %look at this
    file=fullfile(datadir,strcat(sub,'_gordon_parcels_plus_subcort.pconn.nii'));
    subfcmat=cifti_read(file).cdata;
    %eliminate the 7.25s on the diagonal
    for x=1:size(subfcmat)
       subfcmat(x,x)=0;
    end
    %make it 333 X 333
    subfcmat=subfcmat(1:333,1:333);
    %take out the indices that took out for SC-FC sample, for nregions=324
    %subfcmat(indices_to_remove, :) = [];
    %subfcmat(:, indices_to_remove) = [];
    avgweight(n,1)=mean(subfcmat(subfcmat~=0));
    
    %segregation proportionally weighted
    [S, W, B] = segregation(subfcmat,gordon_parcels.numeric_network);
    system_segreg_prop(n,1)=S;
    mean_within_sys_prop(n,1)=W;
    mean_between_sys_prop(n,1)=B;

    %Within and between connectivity, code adapted from Micaela Chan 2018
    Ci=gordon_parcels.numeric_network;
    nCi = unique(Ci);
    M=subfcmat;
    
    for i = 1:length(nCi) % loop through communities
        for j = 1:length(nCi)
           Wi = Ci == nCi(i); % find index for within community edges
           Bi = Ci == nCi(j); % find index for between community edges to specific community
           
           Wv_temp = M(Wi,Wi); % extract within community edges
           Bv_temp = M(Wi,Bi); % extract between community edges to specific community
           
           %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
           Bv = [Bv_temp(:)'];
           system_connectivity(i,j)=mean(Bv(Bv~=0));
           %system_between(i,1)=mean(Bv);
           %if i==j
           %else
        end
    end

%transpose these (or something so they can be saved out on a subject basis.
system_connectivity;
%put this all into a matrix for everyone so that we can see the average
%system connectivity
system_connectivity_all(:,:,n)=system_connectivity;
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

%modularity using Gen Louvain classic modularity maximization, negative
%asymmetric weighting
for c = 1:100 %run it 100x for each subject, get an average number of communities detected.
    [M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
    modul_temp(c)=Q;
    num_comms_temp(c)=length(unique(M));
%also save the number of communities detected.
end
modul(n,1)=mean(modul_temp(:));
num_comms_modul(n,1)=mean(num_comms_temp(:));
%average participation coefficient!
[Ppos, Pneg]=participation_coef_sign(subfcmat, gordon_parcels.numeric_network);
part_coef_pos(n,1)=mean(Ppos);
part_coef_neg(n,1)=mean(Pneg);
%write out all nodes participation coefficient
part_coef_avg_all(n,:)=((Ppos+Pneg)/2);

%average clustering coefficient
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));
%write out all nodes clustering coefficient
avgclustco_all(n,:)=clustering_coef_wu_sign(subfcmat,3);
end

save('~/Box/projects/in_progress/n88_y3_network_metrics.mat')

header={'ID', 'avgweight', 'modul_avg', 'avg_clustco_both','num_comms_modul_avg','part_coef_pos','part_coef_neg', 'system_segreg', 'mean_within_sys', 'mean_between_sys', compose('net%d', 1:78)}

outfile=table(char(subjlist.modid), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg_prop, mean_within_sys_prop, mean_between_sys_prop, system_conn)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, 'n88_y3_within_between_gordon_withmodulpartcoef.mat'), 'outfile2')
writetable(outfile2,fullfile(outdir,'n88_y3_within_between_gordon_withmodulpartcoef.csv'))

%also save the mean system connectivity matrix
mean_system_conn_mat=mean(system_connectivity_all,3)
save(fullfile(outdir, strcat('n88_y3_mean_system_conn.mat')), 'mean_system_conn_mat')

%save the nodal participation coefficient and clustering coefficient
outfile=table(char(subjlist.modid), part_coef_avg_all)
writetable(outfile,strcat(outdir,'/n88_y3_part_coef_avg_nodewise.csv'))

outfile=table(char(subjlist.modid), avgclustco_all)
writetable(outfile,strcat(outdir,'/n88_y3_clust_coef_avg_nodewise.csv'))
