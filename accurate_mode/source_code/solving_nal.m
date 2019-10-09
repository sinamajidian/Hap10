
%fragment_file='data/65a_1_20/frag20.txt';


name_out_mat=frag2mat(fragment_file);

name_hap=strcat(fragment_file(1:length(fragment_file)-3),'hap');

load(name_out_mat)

% tiny 'data/simulation1a/1a_3.mat';
% name_mat=  'data/1m_withN_1a/02/R2.mat';  
% name_mat=  'data/65a_1_20/20.mat';
% 
% name_hap= 'data/1m_withN_1a/02/2_a.hap';%'data/1m_withN_1a/01/nal_hap2_full.hap';%'data/65a_1_20/sdhap_matlab_nal.hap'



addpath(genpath('/Users/sina/Documents/MATLAB/SDPNAL/'))


size(R)
allel_each_row=sum(abs(full(R)),2);    
frags_good_length_ind=allel_each_row>3;
R=R(frags_good_length_ind,:);   
cov=sum(abs(full(R)));
good_cov_snp= find(cov>3); % if coverage of a column is less than 3, remove it.
R=R(:,good_cov_snp);
hap_index=hap_index(good_cov_snp);
allel_each_row=sum(abs(full(R)),2);    
frags_good_length_ind=find(allel_each_row>3);
R=R(frags_good_length_ind,:);  
size(R)


k=3;
N=size(R,1);
l=size(R,2);
W=zeros(N,N);
diag_const=0;
W(1,1)=diag_const;
R1=full(R);
for i=2:N
    line_i=R1(i,:);
    W(i,i)=diag_const;
    for j=1:(i-1)
        line_j=R1(j,:);
        SNP_shared=sum( (line_i~=0) & (line_j~=0));
        allele_shared=sum( line_i.*line_j);
        if SNP_shared>0
            W(i,j)=(2*allele_shared-SNP_shared)/SNP_shared;
        end
    end
end
W=W'+tril(W,-1);
size(W)

 
X=sdp_solver(-W);
 

[Q, sig]=eig(X);
[val_eig, idx]=sort(diag(sig), 'descend'); % ascend[1,2,3]  descend [2,3,1]
three_ind_largest=idx(1:3); % sometimes  the third is zero
val_eig(1:7)'

V=Q(:,three_ind_largest)*sqrt(sig(three_ind_largest,three_ind_largest));


V_arch=V;


object_all=[];
indx_all=[];

 num_it=10*floor(log2(N));
for ii=1:num_it
    Z=normrnd(0,1,[k,k]); %Z=normalize(Z);
    VZ=V*Z;
    [val, index]=max(VZ');
    X_estimated=ones(N,N);
    for i=1:N
        for j=1:N
            if index(i)~=index(j)
                X_estimated(i,j)=-1;
            end
        end
    end
    object_all=[object_all; trace(W*X_estimated)  ];
    indx_all=[indx_all;index];
end

[vall,i_best]=max(object_all);
index_best=indx_all(i_best,:);
 
size(index_best)

R=full(R);
H=zeros(3,l);

for i_k=1:k
   R(index_best==i_k,:);
   H(i_k,:)=sum(R(index_best==i_k,:))>0;
end


% H_final=2*H-1;

% save('sim_R2_nal_full.mat','-v7.3')

% % % % %%% greedy refinement
H_new=2*H-1;


H_final=refiner(R,H_new);





% sum(abs(H_sd(1,:)-H_final(3,:)))
indces_block=hap_index'-1;  % The output file will be like sdhap. index starts from zero
fileID_hap = fopen(name_hap,'w');
fprintf(fileID_hap,'Block 1\t Length of haplotype block %d\t Number of read %d\t Total MEC 01 \n',length(indces_block),N);
H_with_ind=[indces_block, (H_final'+1)/2+1];
fprintf(fileID_hap,'%d\t%d\t%d\t%d\n',H_with_ind');
% 
% % 
% % 


