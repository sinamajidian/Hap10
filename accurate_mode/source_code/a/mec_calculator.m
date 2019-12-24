

function mec=mec_calculator(R_matrix,H_candidate)

% R {+1,-1,0}   h {1,-1}
N=size(R_matrix,1);
[K, ~]=size(H_candidate); % K rows , Kploidy level

mec=0;
for i=1:N
    R_i=R_matrix(i,:);
    [~, R_i_ind,R_i_val]=find(R_i); %[I,J,value]=find() % this is a vector so I=ones
    H_i=H_candidate(:,R_i_ind);
    diff=repmat(R_i_val,K,1)-H_i;
    sm=sum(abs(diff),2)/2;
    mec=mec+min(sm);
end
end
