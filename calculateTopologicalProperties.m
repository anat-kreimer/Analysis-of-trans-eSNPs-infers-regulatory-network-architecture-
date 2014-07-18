PPI_network_pairs=dlmread('PPI_network.txt');
%%% create PPI network in a matrix form %%%
ppi_nodes_list=unique([PPI_network_pairs(:,1)' PPI_network_pairs(:,2)']);
nof_genes_ppi=length(ppi_nodes_list);
PPI_network=zeros(nof_genes_ppi);
for i=1:length(PPI_network_pairs(:,1))
    loc1=(find(ppi_nodes_list==PPI_network_pairs(i,1)));
    loc2=(find(ppi_nodes_list==PPI_network_pairs(i,2)));
    PPI_network(loc1,loc2)=1;
    PPI_network(loc2,loc1)=1;
    clear loc1 loc2
end
%%% Calculate nodes rank and distance between every pair of nodes in the network %%%
INFINITY = 100000;
%Initialize
dist = zeros(nof_genes_ppi);
nodes_rank = zeros(1, nof_genes_ppi);
for i = 1:nof_genes_ppi
    for j = 1:nof_genes_ppi 
        %calculate the rank for each node
        if (PPI_network(i,j)== 1 )
           nodes_rank(i)=nodes_rank(i)+1;
        end               
        if (PPI_network(i,j)== 0)
            if (i ~= j)
                dist(i,j)=INFINITY; 
            end
        else 
            dist(i,j)=1;
            if(i==j)
                dist(i,j)=0;
            end
        end
     end
end
%fprintf('Initialize done!\n');
%Main loop
for k = 1:nof_genes_ppi 
    for i = 1:nof_genes_ppi 
        for j = 1:nof_genes_ppi 
            if (dist(i,j) > dist(i,k)+dist(k,j)) 
                dist(i,j)=dist(i,k)+dist(k,j);
            end
        end
    end
end
second_max=max(max(dist(find(dist<INFINITY))));
for i = 1:nof_genes_ppi 
    for j = 1:nof_genes_ppi 
        if (dist(i,j) == INFINITY)
            dist(i,j) = 2*second_max;	
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Claculate topological properties for assciation pairs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_target_snp_pval=dlmread('source_target_snp_pval.txt');
source_target_dist_sdegree_tdegree=zeros(length(source_target_snp_pval(:,1)),5);
for i=1:length(source_target_snp_pval(:,1))
    source_target_dist_sdegree_tdegree(i,1)=source_target_snp_pval(i,1);
    source_target_dist_sdegree_tdegree(i,2)=source_target_snp_pval(i,2);
    loc_source=(find(ppi_nodes_list==source_target_snp_pval(i,1)));
    loc_target=(find(ppi_nodes_list==source_target_snp_pval(i,2)));
    source_target_dist_sdegree_tdegree(i,3)=dist(loc_source,loc_target);
    source_target_dist_sdegree_tdegree(i,4)=nodes_rank(loc_source);
    source_target_dist_sdegree_tdegree(i,5)=nodes_rank(loc_target);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PRINT %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print topological properties %%%
fid=fopen('topological_properties_of_eSNPs.txt','w+');
for i=1:length(source_target_dist_sdegree_tdegree(:,1))
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',source_target_dist_sdegree_tdegree(i,:));
end
fclose(fid);


