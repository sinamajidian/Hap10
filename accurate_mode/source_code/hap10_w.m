% solving
%clearvars
%i=0;j=2;
%fragment_file=strcat('data/5m_.001_1_cov5/',num2str(i),'/frag',num2str(i),'_',num2str(j),'.txt')

%function H_final=hap10(fragment_file,K)
K=3;
tic
fragment_file= '/mnt/LTR_userdata/majid001/nobackup/5m/4_c10/hap10_2/2/frag2_1.txt'
name_out_mat=frag2mat(fragment_file);

name_hap=strcat(fragment_file(1:length(fragment_file)-3),'hap');
load(name_out_mat)
%cove 10
%min_allele_row=3; 
%min_cov=4;
size(R)
[N, l]=size(R);
cov=sum(abs(full(R)));
if mean(cov)>15 && N>9000
min_allele_row=3;
min_cov=4;

else
min_allele_row=2
min_cov=2
end

%cov2
%min_allele_row=2; 
%min_cov=3;

%min_allele_row=4; 
%min_cov=3;


[mean(cov),  min(cov)];
allel_each_row=sum(abs(full(R)),2);
[mean(allel_each_row), min(allel_each_row)]

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


if (N>2 && l>1)

%W=zeros(N,N);
%diag_const=0;
%W(1,1)=diag_const;
R1=full(R);
%for i=2:N
    %line_i=R1(i,:);
    %W(i,i)=diag_const;
    %for j=1:(i-1)
        %line_j=R1(j,:);
        %SNP_shared=sum( (line_i~=0) & (line_j~=0));
        %#allele_shared=sum( line_i.*line_j); # wrong 
        %allele_shared=sum( (line_i.*line_j)==1); 
        %if SNP_shared>0
            %W(i,j)=(2*allele_shared-SNP_shared)/SNP_shared;
        %end
    %end
%end
%W=W'+tril(W,-1);

toc  # 1
k_sim_k_dis=R1*R1';
R1a=abs(R1);     
SNP_shared_mat=R1a*R1a'; %omeg=find(SNP_shared_mat)
W=k_sim_k_dis./SNP_shared_mat;
W(isnan(W))=0;
toc # 2


X=sdp_solver(-W);

toc #3 after     breakyes ..

[Q, sig]=eig(X);
[val_eig, idx]=sort(diag(sig), 'descend'); % ascend[1,2,3]  descend [2,3,1]
three_ind_largest=idx(1:K); % sometimes  the third is zero
%val_eig(1:4)'

V=Q(:,three_ind_largest)*sqrt(sig(three_ind_largest,three_ind_largest));


%V_arch=V;

toc # 4 
object_all=[];
indx_all=[];
 num_it=50*floor(log2(N));
for ii=1:num_it
    Z=normrnd(0,1,[K,K]); %Z=normalize(Z);
    VZ=V*Z;
    [~, index]=max(VZ'); % no max(VZ,[],2);
    %X_estimated=ones(N,N);
    %for i=1:N
        %for j=1:N
            %if index(i)~=index(j)
                %X_estimated(i,j)=-1;
            %end
        %end
    %end
    index_mat=repmat(index,N,1);
    X_estimated=2*(index_mat==index_mat')-1;
    object_all=[object_all; trace(W*X_estimated)  ];
    indx_all=[indx_all;index];
end
toc
[~,i_best]=max(object_all);
index_best=indx_all(i_best,:);


R=full(R);
H=zeros(K,l);

for i_k=1:K
   R(index_best==i_k,:);
   H(i_k,:)=sum(R(index_best==i_k,:))>0;
end



% % % % %%% greedy refinement
H_new=2*H-1;

toc
H_final=refiner(R,H_new);
toc
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
toc

%end
