% solving
%clearvars



%fragment_file='/mnt/LTR_userdata/majid001/nobackup/1m/1m01/1_c5/hap10_2/0/frag0_1.txt'
%K=3;


function H_final=hap10(fragment_file,K)

tic
name_out_mat=frag2mat(fragment_file);

name_hap=strcat(fragment_file(1:length(fragment_file)-3),'hap');
load(name_out_mat)


cov=sum(abs(full(R)));
size(R)
[N, l]=size(R);

if mean(cov)>15 && N>9000
min_allele_row=3;
min_cov=3;

else
min_allele_row=2;
min_cov=0;
end


[mean(cov),  min(cov)];
allel_each_row=sum(abs(full(R)),2);
[mean(allel_each_row), min(allel_each_row)];

allel_each_row=sum(abs(full(R)),2);    
frags_good_length_ind=allel_each_row>=min_allele_row;
R=R(frags_good_length_ind,:);   
cov=sum(abs(full(R)));
good_cov_snp= find(cov>=min_cov); % if coverage of a column is less than 3, remove it.
R=R(:,good_cov_snp);
hap_index=hap_index(good_cov_snp);
allel_each_row=sum(abs(full(R)),2);    
frags_good_length_ind=find(allel_each_row>=min_allele_row);
R=R(frags_good_length_ind,:);  

%cov=sum(abs(full(R)));
%[mean(cov),  min(cov)];
%allel_each_row=sum(abs(full(R)),2);
%[mean(allel_each_row), min(allel_each_row)];


size_R=size(R)
N=size(R,1);
l=size(R,2);


if (N>2 && l>3)

%W=zeros(N,N);
%diag_const=0;
%W(1,1)=diag_const;
R1=full(R);

k_sim_k_dis=R1*R1';
R1a=abs(R1);     
SNP_shared_mat=R1a*R1a'; %omeg=find(SNP_shared_mat)
W=k_sim_k_dis./SNP_shared_mat;
W(isnan(W))=0;


toc %1
X=sdp_solver(-W);
toc  %2
[Q, sig]=eig(X);
[val_eig, idx]=sort(diag(sig), 'descend'); % ascend[1,2,3]  descend [2,3,1]
three_ind_largest=idx(1:K); % sometimes  the third is zero
%val_eig(1:4)'

V=Q(:,three_ind_largest)*sqrt(sig(three_ind_largest,three_ind_largest));


%V_arch=V;

toc %3


num_it=1000*floor(log10(N)); % 100
object_all=zeros(num_it,1);
indx_all=zeros(num_it,N);
    W_sp=sparse(W);
for ii=1:num_it
    Z=normrnd(0,1,[K,K]); 
    VZ=V*Z;
    [~, index]=max(VZ'); % no max(VZ,[],2);
    index_mat=repmat(index,N,1);
    X_estimated=2*(index_mat==index_mat')-1;    
    object_all(ii)= W_sp(:).'*reshape(X_estimated.',[],1); %trace(W_sp*X_estimated);
    indx_all(ii,:)=index;
end
[~,i_best]=max(object_all);
index_best=indx_all(i_best,:);


R=full(R);
H_b=zeros(K,l);
HQ=zeros(K,l);


for i_k=1:K
  % R(index_best==i_k,:);
   value=sum(R(index_best==i_k,:));
   H_b(i_k,:)=value>0;
   HQ(i_k,:)=abs(value);
end
H_one=2*H_b-1;



% % % % %%% greedy refinement
H_new=H_one;

toc %4
H_final=refiner(R,H_new);
toc %5
mec_final=mec_calculator(R,H_final);




indces_block=hap_index'-1;  % The output file will be like sdhap. index starts from zero
H_with_ind=[indces_block, (H_final'+1)/2+1];

fileID_hap = fopen(name_hap,'w'); 
fprintf(fileID_hap,'Block 1\t Length of haplotype block %d\t Number of read %d\t Total MEC %d \n',length(indces_block),N,mec_final);
string_d=strcat(repmat('%d\t', 1, K),'%d\n');
fprintf(fileID_hap,string_d,H_with_ind');
fclose(fileID_hap);



name_out_all=strcat(fragment_file(1:length(fragment_file)-4),'_all.mat');
save(name_out_all,'-v7.3')
else
    fprintf(' The number of reads or SNPs is not enough. At least three (proper) reads and two SNPs.  \n')
    H_final=0;
end


end
