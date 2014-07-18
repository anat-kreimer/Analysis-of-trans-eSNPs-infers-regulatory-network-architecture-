%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Find the significant tripltes %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero
Expression_matrix=dlmread('Expression_matrix.txt');
nof_ind=size(Expression_matrix,1);
nof_genes=size(Expression_matrix,2);
SNPs_linear=dlmread('SNPs_minor_allele_count.txt');
nof_snps=size(SNPs_linear,2);
tmp_gene_snp_pval=dlmread('gene_snp_pval.txt');
tmp_gene_snp_pval=unique(tmp_gene_snp_pval,'rows');
[B,I]=sort(tmp_gene_snp_pval(:,2));
gene_snp_pval=tmp_gene_snp_pval(I,:);
clear B I
current_three_some=1;
%%%%%% create SNP_count %%%%%%%%%%
SNPs_count=zeros(nof_snps,1);
for i=1:nof_snps
    if(find(gene_snp_pval(:,2)==i))
       SNPs_count(i)=length(find(gene_snp_pval(:,2)==i));
    end
end
%triplets matrix (#triplets X 4): Gene1, Gene2, SNP
nof_triplets=0;
for i=1:nof_snps
    if(find(gene_snp_pval(:,2)==i))
        if (SNPs_count(i) > 1)
            if(length(find(SNPs_linear(:,i)==0))< nof_ind)
            for x=1:(SNPs_count(i)-1)
                for y =(x+1):SNPs_count(i)
                    %side 1
                    Model_gene1_snp=[Expression_matrix(:,gene_snp_pval(current_three_some+y-1,1)) SNPs_linear(:,i)];
                    [b,dev,stat1] = glmfit(Model_gene1_snp,Expression_matrix(:,gene_snp_pval(current_three_some+x-1,1)));
                    %side 1
                    Model_gene2_snp=[Expression_matrix(:,gene_snp_pval(current_three_some+x-1,1)) SNPs_linear(:,i)];
                    [b,dev,stat2] = glmfit(Model_gene2_snp,Expression_matrix(:,gene_snp_pval(current_three_some+y-1,1)));
                    % Triplet topology
                    if (stat1.p(3)<=0.05 && stat2.p(3)<=0.05)
                        nof_triplets = nof_triplets+1;
                        triplets(nof_triplets,1)=gene_snp_pval(current_three_some+x-1,1);
                        triplets(nof_triplets,2)=gene_snp_pval(current_three_some+y-1,1);
                        triplets(nof_triplets,3)=i;
                    else
                        if (stat1.p(3)<=0.05 && stat2.p(3)>0.05)
                            nof_triplets = nof_triplets+1;
                            triplets(nof_triplets,1)=gene_snp_pval(current_three_some+x-1,1);
                            triplets(nof_triplets,2)=gene_snp_pval(current_three_some+y-1,1);
                            triplets(nof_triplets,3)=i;
                        end
                        if (stat1.p(3)>0.05 && stat2.p(3)<=0.05)
                            nof_triplets = nof_triplets+1;
                            triplets(nof_triplets,1)=gene_snp_pval(current_three_some+y-1,1);
                            triplets(nof_triplets,2)=gene_snp_pval(current_three_some+x-1,1);
                            triplets(nof_triplets,3)=i;
                        end
                    end
                    clear stat1 stat2 b dev
                end
            end
            current_three_some = current_three_some+SNPs_count(i);
            else
                current_three_some=current_three_some+SNPs_count(i);
            end 
        else
            current_three_some=current_three_some+1;
        end
    end
end
clear x y i
for i=1:length(triplets(:,1))
    if(triplets(i,1)>triplets(i,2))
        tmp_var=triplets(i,1);
        triplets(i,1)=triplets(i,2);
        triplets(i,2)=tmp_var;
        clear tmp_var
    end
end
triplets=unique(triplets,'rows');
[B,tmp_ind] = sort(triplets(:,1));
triplets=triplets(tmp_ind,:);
clear B tmp_ind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Filter and Assemble cooperating quartets %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unique_gene_pairs]=unique(triplets(:,[1 2]),'rows');
%go over each pair of SNPs
count=0;
first_filter_2genes_2snps=[];
for i=1:length(unique_gene_pairs(:,1))
    tmp_loc=find(triplets(:,1)==unique_gene_pairs(i,1) & triplets(:,2)==unique_gene_pairs(i,2));
    if(length(tmp_loc)>1)
        for x=1:length(tmp_loc)-1
            for y=(x+1):length(tmp_loc)
                Model_snps=[SNPs_linear(:,triplets(tmp_loc(x),3)) SNPs_linear(:,triplets(tmp_loc(y),3))];
                [b,dev,st1]=glmfit(Model_snps,Expression_matrix(:,triplets(tmp_loc(x),1)));
                [b,dev,st2]=glmfit(Model_snps,Expression_matrix(:,triplets(tmp_loc(x),2)));
                %First filter: cooperating SNPs - all betas are significant
                if(st1.p(2)<=0.05 && st1.p(3)<=0.05 && st2.p(2)<=0.05 && st2.p(3)<=0.05)
                    count=count+1;
                    first_filter_2genes_2snps(count,:)=[triplets(tmp_loc(x),1) triplets(tmp_loc(x),2) triplets(tmp_loc(x),3) triplets(tmp_loc(y),3)];
                end
                clear Model_snps b dev st1 st2
            end
        end
    end
    clear tmp_loc
end
if(first_filter_2genes_2snps)
    [x,y,z]=unique(first_filter_2genes_2snps(:,[1 2]),'rows');
    first_filter_2genes_2snps=first_filter_2genes_2snps(y,:);
else
    fprintf('No quartets found\n');
    exit
end
clear x y z i count

snps_data=dlmread('SNPs_data.txt');
list_of_potential_intermidiate_SNPs=dlmread('potential_intermidiate_SNPs.txt');
%Second filter: no third SNP can explain the expression better
count=0;
excluded_quartets_index=[];
for i=1:length(first_filter_2genes_2snps(:,1))
    % cooperating SNPs are on the same chr
    if (snps_data(first_filter_2genes_2snps(i,3),1)==snps_data(first_filter_2genes_2snps(i,4),1))
        for j=1:length(list_of_potential_intermidiate_SNPs)
            %is the potential SNP on the same chr
            if(first_filter_2genes_2snps(i,3)==snps_data(list_of_potential_intermidiate_SNPs(j),1))
                [r1,t,p]=spear(SNPs_linear(:,first_filter_2genes_2snps(i,3)),SNPs_linear(:,list_of_potential_intermidiate_SNPs(j)));
                [r2,t,p]=spear(SNPs_linear(:,first_filter_2genes_2snps(i,4)),SNPs_linear(:,list_of_potential_intermidiate_SNPs(j)));
                %both r1^2 and r2^2 should be >= 0.5
                if((r1^2 >=0.5) && (r2^2 >=0.5))
                    %check if this intermiditae SNP is in a triplet with the
                    %two genes
                    Model_gene1_snp=[Expression_matrix(:,first_filter_2genes_2snps(i,1)) SNPs_linear(:,list_of_potential_intermidiate_SNPs(j))];
                    Model_gene2_snp=[Expression_matrix(:,first_filter_2genes_2snps(i,2)) SNPs_linear(:,list_of_potential_intermidiate_SNPs(j))];
                    [b,dev,st1]=glmfit(Model_gene1_snp,Expression_matrix(:,first_filter_2genes_2snps(i,2)));
                    [b,dev,st2]=glmfit(Model_gene2_snp,Expression_matrix(:,first_filter_2genes_2snps(i,1)));                   
                    %If one of the betas for the intermidiate SNP is
                    %significant, then the SNP is in a triplet
                    if(st1.p(3)<=0.05 || st2.p(3)<=0.05)
                         Model_3SNPs=[first_filter_2genes_2snps(i,3) first_filter_2genes_2snps(i,4) SNPs_linear(:,list_of_potential_intermidiate_SNPs(j))];
                         [b,dev,st_gene1]=glmfit(Model_3SNPs,Expression_matrix(:,first_filter_2genes_2snps(i,1)));
                         [b,dev,st_gene2]=glmfit(Model_3SNPs,Expression_matrix(:,first_filter_2genes_2snps(i,2)));
                         % does SNP3 explain better???
                         if(st_gene1.p(4) <= 0.05 && st_gene2.p(4) <= 0.05)
                             count=count+1;
                             excluded_quartets_index(count)=i;
                         end
                         clear Model_3SNPs st_gene1 st_gene2 b dev
                    end
                    clear Model_gene1_snp Model_gene2_snp st1 st2 b dev            
                end
            end
            clear r1 r2 t p 
        end
    end
end
clear i j
count=0;
for i=1:length(first_filter_2genes_2snps(:,1))
    if(excluded_quartets_index)
        if(find(excluded_quartets_index==i))
        else
            count=count+1;
            cooperating_quartets(count,:)=[first_filter_2genes_2snps(i,1) first_filter_2genes_2snps(i,2) first_filter_2genes_2snps(i,3) first_filter_2genes_2snps(i,4)];
        end
    else
        cooperating_quartets(i,:)=[first_filter_2genes_2snps(i,1) first_filter_2genes_2snps(i,2) first_filter_2genes_2snps(i,3) first_filter_2genes_2snps(i,4)];
    end
end
clear count
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PRINT %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print cooperating quartets %%%
fid=fopen('quartets.txt','w+');
for i=1:length(cooperating_quartets(:,1))
    fprintf(fid,'%d\t%d\t%d\t%d\n',cooperating_quartets(i,:));
end
fclose(fid);

   
    
