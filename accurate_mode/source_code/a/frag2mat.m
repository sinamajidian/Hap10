
function name_out=frag2mat(fragment_file)

name_out=strcat(fragment_file(1:length(fragment_file)-3),'mat');

% please check that framgent starts from first row or third row



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Converting a fragment file (ProbHAP) to sparse matrix


% Input: a fragment file
% output: a .mat file containg sparse matrix


% a sample fragment file
% 2         % Number of reads  % we dont use this line, the first line  is discarded by readtable() 
% 200       % Number of columns 
% 1 NC_001133.9_1-318 166 0000101001 AIGGGGDIII
% 1 NC_001133.9_2-288 166 0000010000 ?HHHEGFGHI

% a sample .mat file
% [-1,-1,1,0,0,0,0,0;
%  -1,-1,-1,-1,0,0,0,0;
%  -1,-1,-1,-1,0,0,0,0]

%Sina Majidian Dec 2018
%Iran University of Science and Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%fid = fopen(fragment_file);
%file = textscan(fid,'%s','delimiter','\n');

%full_table=readtable(fragment_file,'Delimiter','\t'); %the first line is ignored
full_table=readtable(fragment_file,'Delimiter','\t','ReadVariableNames',false);
full_cell=table2cell(full_table); 
%l=str2double(full_cell{2});       % number of SNPs i.e. haplotype length
fragment_cell =full_cell;%(1:end); % removing two first line (number of read and number of snps=haplotype length )

[first_index,last_index]=calculate_first_last(fragment_cell);
hap_index=first_index:last_index;
l=last_index; % this l is not haplotype length.  first index may not be  zero 
start_i=first_index;

N=size(fragment_cell,1);   % the number of reads i.e. row in fragment file
R=sparse(N,l); % the final read matrix


for i=1:N % index for each row of the cell (file)
    row_num=zeros(1,l);
    
    %1 chr22_SPH2_1940 3 100 C==
    % string is an element of cell,
    %row_str = fragment_cell{N-i+1};  %% for haplogenerator
    row_str = fragment_cell{i};      %% for dutima data, sorted from first row
    
    
    
    ii=1; % index of each letter in a row
    % each row may contain several blocks. Frist number in a row is the number of blocks,
    String_BlockNumber=[];        %number of blocks as a a string format
    while isspace(row_str(ii))~=1 %number of blocks, it can be mor than one digit
        String_BlockNumber=[String_BlockNumber,row_str(ii)];
        ii=ii+1;
    end
    BlockNumber=str2double(String_BlockNumber);
    startingPoint_num=zeros(1,BlockNumber); %A vector of starting points. For each block, a starting point exist
    ii=ii+1;
    while isspace(row_str(ii))~=1 % for going after read id 'Chr2_...'
        ii=ii+1;
    end
    for j=1:BlockNumber
        ii=ii+1;
        startingPoint_str=[]; %each starting point as a string
        while isspace(row_str(ii))~=1
            startingPoint_str=[startingPoint_str,row_str(ii)];
            ii=ii+1;
        end
        startingPoint_num(j)=str2double(startingPoint_str);
        site=startingPoint_num(j); % index of each row of final matrix
        ii=ii+1;
        while isspace(row_str(ii))~=1  %
            row_num(site)=2*str2double(row_str(ii))-1;%put the values {0,1} in 'row'  for each block as {1,-1}
            site=site+1;
            ii=ii+1;
        end
    end
    R(i,:)=row_num; % insert this row as the row of final matrix
end


R=R(:,start_i:end);
R1=full(R);



% clearvars -except R  fragment_cell name_out hap_index
save(name_out,'-v7.3')

%R_f=full(R);








function [first_index_global,last_index_global]=calculate_first_last(fragment_cell)

N=size(fragment_cell,1);   % the number of reads i.e. row in fragment file
first_index_global=1e9;
last_index_global=0;

 for i=1:N % index for each row of the cell (file)
    row_str = fragment_cell{i};      %% for dutima data, sorted from first row  
    ii=1; % index of each letter in a row
    % each row may contain several blocks. Frist number in a row is the number of blocks,
    String_BlockNumber=[];   %number of blocks as a a string format
    while isspace(row_str(ii))~=1 %number of blocks, it can be mor than one digit
        String_BlockNumber=[String_BlockNumber,row_str(ii)];
        ii=ii+1;
    end
    BlockNumber=str2double(String_BlockNumber);
    startingPoint_num=zeros(1,BlockNumber); %A vector of starting points. For each block, a starting point exist
    ii=ii+1;
    while isspace(row_str(ii))~=1 % for going after read id 'Chr2_...'
        ii=ii+1;
    end
    
    for j=1:BlockNumber
        ii=ii+1;
        startingPoint_str=[]; %each starting point as a string
        while isspace(row_str(ii))~=1
            startingPoint_str=[startingPoint_str,row_str(ii)];
            ii=ii+1;
        end
        startingPoint_num(j)=str2double(startingPoint_str);
        site=startingPoint_num(j); % index of each row of final matrix
        ii=ii+1;
        while isspace(row_str(ii))~=1  %
            row_num(site)=2*str2double(row_str(ii))-1;%put the values {0,1} in 'row'  for each block as {1,-1}
            site=site+1;
            ii=ii+1;
        end
        last_index_global= max( [last_index_global ,site-1]);
    end
    first_index_global=min([first_index_global,startingPoint_num ]);
 end
end



end
