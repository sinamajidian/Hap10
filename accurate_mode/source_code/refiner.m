


function H_final=refiner(R,H)
% element of the input H should be 1.-1

[K, L]=size(H); % k rows , k ploidy level

H_best_upnow=H;
KK=10;
mec_all=zeros(1,KK);

mec_all(1)=mec_calculator(R,H_best_upnow);
kk=1;
continu=1;
while (kk<KK && continu)  % I can save all of them and extract the best mec, no it can become worde
    kk=kk+1;
    
    for k=1:K % row  of H
        for l=1:L % col of H     
            H_check=H_best_upnow;
            H_check(k,l) = -H_best_upnow(k,l);
            mec_check_l=mec_calculator_l(R,H_check,l);
            mec_best_l=mec_calculator_l(R,H_best_upnow,l);
            if mec_check_l < mec_best_l
                H_best_upnow=H_check;
            end
        end
    end


mec_all(kk) = mec_calculator(R,H_best_upnow);

if ( (mec_all(kk)==mec_all(kk-1)) || mec_all(kk)==0 )
    continu=0;
end

mec_all_end=mec_all(1:kk);

H_final=H_best_upnow;

end

mec_all_end







function mec=mec_calculator_l(R_matrix,H_candidate,l)

% column l in read matrix, l-th snp position in haplotyp.
% R {+1,-1,0}   h {1,-1}
%R_matrix_l=R_matrix(:,l);
[K, ~]=size(H_candidate); % K rows , K ploidy level

[row_list,~,~]=find(R_matrix(:,l)); % row list that are non zero in l-th column
mec=0;
%R_rows_inlist=R_matrix(row_list,:);
for i=1:length(row_list)
    %R_i=R_rows_inlist(i,:);
    R_i=R_matrix(row_list(i),:);

    [~, R_i_ind, R_i_val]=find(R_i); %[I,J,value]=find() % this is a vector so I=ones
    H_i=H_candidate(:,R_i_ind);
    diff=repmat(R_i_val,K,1)-H_i;
    sm=sum(abs(diff),2)/2;
    mec=mec+min(sm);
end
end


end
