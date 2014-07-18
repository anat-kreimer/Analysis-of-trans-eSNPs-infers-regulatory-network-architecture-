%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Find the significant tripltes %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero
Expression_matrix=dlmread('Expression_matrix.txt');
nof_ind=size(Expression_matrix,1);
nof_genes=size(Expression_matrix,2);
SNPs_linear=dlmread('SNPs_minor_allele_count.txt',' ');
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
%triplets matrix (#triplets X 4): Gene1, Gene2, SNP, pval1, pval2, edge type (Gene1--->Gene2; 1 - uni, 2 -bi)
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
                        triplets(nof_triplets,4)=-log10(gene_snp_pval(current_three_some+x-1,3));
                        triplets(nof_triplets,5)=-log10(gene_snp_pval(current_three_some+y-1,3));
                        triplets(nof_triplets,6)=2;
                    else
                        if (stat1.p(3)<=0.05 && stat2.p(3)>0.05)
                            nof_triplets = nof_triplets+1;
                            triplets(nof_triplets,1)=gene_snp_pval(current_three_some+x-1,1);
                            triplets(nof_triplets,2)=gene_snp_pval(current_three_some+y-1,1);
                            triplets(nof_triplets,3)=i;
                            triplets(nof_triplets,4)=-log10(gene_snp_pval(current_three_some+x-1,3));
                            triplets(nof_triplets,5)=-log10(gene_snp_pval(current_three_some+y-1,3));
                            triplets(nof_triplets,6)=1;
                        end
                        if (stat1.p(3)>0.05 && stat2.p(3)<=0.05)
                            nof_triplets = nof_triplets+1;
                            triplets(nof_triplets,1)=gene_snp_pval(current_three_some+y-1,1);
                            triplets(nof_triplets,2)=gene_snp_pval(current_three_some+x-1,1);
                            triplets(nof_triplets,3)=i;
                            triplets(nof_triplets,4)=-log10(gene_snp_pval(current_three_some+y-1,3));
                            triplets(nof_triplets,5)=-log10(gene_snp_pval(current_three_some+x-1,3));
                            triplets(nof_triplets,6)=1;
                        end
                    end
                    clear stat1 stat2 b dev Model_gene1_snp Model_gene2_snp
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
    fprintf('%d\n',i);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Assemble modules %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snps_modules=unique(triplets(:,3));
nof_modules=length(snps_modules);
snp_index_module_size=zeros(nof_modules,2);
for i=1:nof_modules
    tmp_vec = find(triplets(:,3)== snps_modules(i));
    len = length(tmp_vec);
    snp_index_module_size(i,1)= snps_modules(i);
    snp_index_module_size(i,2)= length(unique(triplets(tmp_vec(1):tmp_vec(len),1:2)));
    clear tmp_vec len
end
%sort by modules size (#genes)
[c,index]=sort(snp_index_module_size(:,2));
snp_index_module_size_sorted=zeros(nof_modules,2);
for i=1:nof_modules
    snp_index_module_size_sorted(i,1)=snp_index_module_size(index(i),1);
    snp_index_module_size_sorted(i,2)=snp_index_module_size(index(i),2);
    tmp_vec = find(triplets(:,3)== snp_index_module_size_sorted(i,1));
    len = length(tmp_vec);
    modules_sorted{i}=unique(triplets(tmp_vec(1):tmp_vec(len),1:2));
    pvals_modules_sorted{i}=cat(1,triplets(tmp_vec(1):tmp_vec(len),4), triplets(tmp_vec(1):tmp_vec(len),5));
    clear tmp_vec len
end
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PRINT %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print main SNP index %%%
fid=fopen('main_SNP.txt','w+');
for i=1:nof_modules
    fprintf(fid,'%d\n',snp_index_module_size_sorted(i,1));
end
fclose(fid);
%%% print genes indices for each module %%%
fid=fopen('modules_sorted.txt','w+');
for i=1:nof_modules
    for j=1:length(modules_sorted{i})
        fprintf(fid,'%d\t',modules_sorted{i}(j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%% print modules sizes %%%
fid=fopen('modules_sizes.txt','w+');
for i=1:nof_modules
    fprintf(fid,'%d\n',snp_index_module_size_sorted(i,2));
end
fclose(fid);
%%% print pvals for each module %%%
fid=fopen('pvals_modules.txt','w+');
for i=1:nof_modules
    for j=1:length(pvals_modules_sorted{i})
        fprintf(fid,'%f\t',pvals_modules_sorted{i}(j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
clear i j fid

   
    
