

K=4; 


fragment_files={};
for cov=[10,2]

for i1=1:5  


K=4; 

adress_folder=strcat('/mnt/LTR_userdata/majid001/nobackup/1m/tet_1m01/',num2str(i1),'_c',num2str(cov),'/hap10_2/');  

all_files = dir(adress_folder);
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir)-2;



for i=0:num_dir-1
all_files_i = dir(strcat(adress_folder,num2str(i)));
file_num=length({all_files_i.name});
frags_num=file_num; 
for j=1:frags_num

fragment_file=strcat(adress_folder,num2str(i),'/frag',num2str(i),'_',num2str(j),'.txt');

fragment_files{end+1}=fragment_file;

end
end

end

%end



size(fragment_files,2) 
cpu_time_all=0;
for i=1:size(fragment_files,2)
fragment_file=fragment_files{i};
if isfile(fragment_file)
    fragment_file
    start_time=cputime;
    H_final=hap10(fragment_file,K);
    end_time=cputime;
    cpu_time=end_time-start_time
    cpu_time_all=cpu_time_all+cpu_time;
end

end

cpu_time_all


#end


%I=[0 0 0 1 1 2 2 2]
%J=[6 7 8 3 4 7 8 9]


%for k=1:length(I)
%i=I(k);
%j=J(k);

%i=3;
%j=1;

%fragment_file=strcat('data/5m_.01_cov5_rerun/',num2str(i),'/frag',num2str(i),'_',num2str(j),'.txt')
%H_final=hap10(fragment_file,K);
%%end
