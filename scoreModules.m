%%%%% set #permutations %%%%%%
nof_perm=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assign a score for each     %%%
%%% module in the permutations  %%%
%%% using a "leave one out"     %%%
%%% procedure                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modules_scores_perm=cell(nof_perm,1);
modules_pvals_fp=cell(nof_perm,1);
modules_fp_rate=cell(nof_perm,1);
%read the permutations modules data into variables%
count_pvals=1;
for x=1:nof_perm
    file_name=sprintf('associations_perm_%d.txt',x);
    tmp=dlmread(file_name);
    perm_all_pvals(count_pvals:1:count_pvals+length(tmp(:,3))-1,1)=-log10(tmp(:,3));
    perm_all_pvals(count_pvals:1:count_pvals+length(tmp(:,3))-1,2)=x;
    count_pvals=count_pvals+length(tmp(:,3));
    file_name=sprintf('modules_sizes_perm_%d.txt',x);
    permutations_modules_sizes{x}=dlmread(file_name);
    file_name=sprintf('pvals_modules_perm_%d.txt',x);
    permutations_modules_pvals{x}=dlmread(file_name);
    clear tmp
end    
for x=1:nof_perm
    fprintf('nof_perm=%d\n',x);
    tmp_nof_modules=size(permutations_modules_sizes{x},2);
    modules_pvals_fp{x}=cell(tmp_nof_modules,1);
    modules_fp_rate{x}=zeros(tmp_nof_modules,1);
    modules_scores_perm{x}=zeros(tmp_nof_modules,1);
    tmp_leave_one_out_loc=find(perm_all_pvals(:,2)~=x);
    tmp_all_but_one_pvals=perm_all_pvals(tmp_leave_one_out_loc,1);
    tmp_one_perm_loc=find(perm_all_pvals(:,2)==x);
    tmp_one_perm_pvals=perm_all_pvals(tmp_one_perm_loc,1);
    for i=1:tmp_nof_modules
        y=1;
        tmp_size=permutations_modules_sizes{x}(i);
        for j=1:tmp_size
            m=length(find(tmp_all_but_one_pvals>=permutations_modules_pvals{x}(i,j)));
            m=m/(nof_perm-1);
            n=length(find(tmp_one_perm_pvals>=permutations_modules_pvals{x}(i,j)));
            modules_pvals_fp{x}{i}(j)= m/n;
            if (n==0)
                modules_pvals_fp{x}{i}(j)= 1;
            end
            if (modules_pvals_fp{x}{i}(j)>1)
                modules_pvals_fp{x}{i}(j)= 1;
            end
            y=y*modules_pvals_fp{x}{i}(j);
        end
        s=0;
        for z=1:nof_perm
            if (z~=x)
                s=s+length(find(permutations_modules_sizes{z}>=permutations_modules_sizes{x}(i)));
            end
        end
        s=s/(nof_perm-1);
        t=length(find(permutations_modules_sizes{x}>=permutations_modules_sizes{x}(i)));
        modules_fp_rate{x}(i)= s/t;
        if (t==0)
            modules_fp_rate{x}(i)=1;
        end   
        if (modules_fp_rate{x}(i)>1)
           modules_fp_rate{x}(i)= 1;
        end
        %module's score
        modules_scores_perm{x}(i)=-log10(y*modules_fp_rate{x}(i));
    end
    clear perm_leave_one_out_pvals tmp_nof_modules tmp_leave_one_out_loc tmp_all_but_one_pvals tmp_one_perm_loc tmp_one_perm_pvals
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Score modules %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modules_pvals_fp=cell(nof_modules,1);
modules_fp_rate=zeros(nof_modules,1);
modules_scores=zeros(nof_modules,1);
for i=1:nof_modules
    fprintf('module number %d\n',i);
    tmp_score=1;
    for j=1:length(modules_sorted{i})
        x=length(find(perm_all_pvals(:,1)>= pvals_modules_sorted{i}(j)));
        x=x/nof_perm;
        y=length(find(-log10(gene_snp_pval(:,3))>= pvals_modules_sorted{i}(j)));
        modules_pvals_fp{i}(j)= x/y;
        tmp_score=tmp_score*modules_pvals_fp{i}(j);
    end
    s=0;
    for z=1:nof_perm
        s=s+length(find(permutations_modules_sizes{z}>=snp_index_module_size_sorted(i,2)));
    end
        s=s/nof_perm;
        t=length(find(snp_index_module_size_sorted(:,2)>=snp_index_module_size_sorted(i,2)));
        modules_fp_rate(i)= s/t;
        %Module's score
        modules_scores(i)=-log10(tmp_score*modules_fp_rate(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PRINT %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print modules' score %%%
fid=fopen('modules_scores.txt','w+');
for i=1:nof_modules
    fprintf(fid,'%f\n',modules_scores(i));
end
fclose(fid);
    



    
        
