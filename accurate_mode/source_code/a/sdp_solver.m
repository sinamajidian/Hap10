
function X_out=sdp_solver(W)


addpath(genpath('/mnt/LTR_userdata/majid001/software/code_matlab/SDPNAL'))

N=size(W,1);

blk{1,1} = 's'; blk{1,2} = N;
C=cell(1);
L=cell(1);
At=cell(1);

C{1}=sparse(W);

L{1}=-.5*ones(N); % L<= X <= U 
U{1}=ones(N);


A_all=sparse((N+1)*N/2,N);
for i=1:N
    ind_i=i*(i+1)/2;
    A_all(ind_i,i)=1;
end

At{1}= A_all; 

b=ones(1,N);

%X_init=cell(1);
%X_init{1}= ones(N);


OPTIONS.printlevel=0; % 1 half report, 2 full report
OPTIONS.tol=1; % 1e-2
OPTIONS.ADMtol=1;
OPTIONS.AATsolve.method='iterative'; % 'direct' (default) 
OPTIONS.BBTsolve.method='iterative';

OPTIONS.ADMmaxiter=400; 
OPTIONS.maxiter=400;
OPTIONS.maxtime=1000; % 10000 (default)  for each of ADM (init) and solve

%[obj,X,s,y,Z1,Z2,y2,v,info,runhist]=sdpnalplus(blk,At,C,b,L,[],[],[],[],OPTIONS,X_init);

[obj,X,s,y,Z1,Z2,y2,v,info,runhist]=sdpnalplus(blk,At,C,b,L,[],[],[],[],OPTIONS,[]);


X_out=X{1};
end
