#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <atlas_enum.h>
#include <clapack.h>



/* DGESDD prototype */
extern void sgesdd_(char* jobz, int* m, int* n, float* a,
int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt,
float* work, int* lwork, int* iwork, int* info);
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, float* a, int lda );
#define Max_row 1000000
#define Max_col 1000000



int main(int argc, char *argv[])
{
	if (argc<2){ // How to use this code :  ./exec_file frag_file output_hap number_haps
		printf("Error: The number of the arguments is %d which should be 3. \n",argc-1 );
		exit(-1);
	}
	FILE *fpt,*fpt1,*fpt2,*fpt5,*fpt7, *file_cc_pt;
	int *row;
	int *col;
	int *val;
	int i,i2;
	int max_row, max_col;
	int iter,iter2;
	int cnt;
	clock_t t1;

	char tempo[100000];
	char qlty[100000];
	char str[100000];
	char *ru;
	int num_read;
	int num_col;
	int **snp_mat;
	int **snp_list;
	int *snp_ct;
	int MEC[100000];
	int num_base[100000];
	float SWER[100000];
	int block_len[100000];
	int block_rdlen[100000];
	int no_acc=0;

	int **arr_pos;
	arr_pos=malloc(100000*sizeof(int*));
	for (iter=0;iter<100000;iter++){
		arr_pos[iter]=malloc(100000*sizeof(int));
	}

	t1=clock();

	int K;
	fpt=fopen(argv[1],"r");
	fpt1=fopen(argv[1],"r");
	//K=atoi(argv[3]);
	fscanf(fpt,"%d",&num_read);
	fscanf(fpt,"%d",&num_col);
	fscanf(fpt1,"%d",&num_read);
	fscanf(fpt1,"%d",&num_col);
	snp_mat=malloc(num_read*sizeof(int*));
	snp_list=malloc(num_read*sizeof(int*));
	snp_ct=malloc(num_read*sizeof(int));
	int **true_hap;
	int temp;
	true_hap=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		true_hap[iter]=malloc(num_col*sizeof(int));
	}
	int ct;
	int num_seg;
	int k;
	int pos;
	int d;
	int j;
	int *col_ct;
	int col_ct2;
	int *rd_array;
	int *rd_array2;
	int **list, **list2;
	int **list_val, **list_val2;
	int *num_ct, *num_ct2;
	int *read_ac;
	i=-1;i2=-1;
	num_ct=malloc(Max_row*sizeof(int));
	read_ac=malloc(Max_row*sizeof(int));
	list=malloc(Max_row*sizeof(int*));
	list_val=malloc(Max_row*sizeof(int*));
	for (iter=0;iter<Max_row;iter++){
		list[iter]=malloc(10000*sizeof(int));
		list_val[iter]=malloc(10000*sizeof(int));
	}
	num_ct2=malloc(Max_col*sizeof(int));
	list2=malloc(Max_col*sizeof(int*));
	list_val2=malloc(Max_col*sizeof(int*));
	for (iter=0;iter<Max_col;iter++){
		list2[iter]=malloc(10000*sizeof(int));
		list_val2[iter]=malloc(10000*sizeof(int));
	}
	rd_array=malloc(num_read*sizeof(int));
	rd_array2=malloc(num_read*sizeof(int));
	col_ct2=0;
	col_ct=malloc(num_col*sizeof(int));
	for(k=0;k<num_col;k++){
		col_ct[k]=0;
	}
	if(fpt==NULL)
	printf("\nError\n");
	else {
		while(i<num_read-1){ // while(!feof(fpt)){
			ct=-1;
			i++;
			fscanf(fpt,"%d",&num_seg);
			fscanf(fpt,"%s",tempo);
			fscanf(fpt1,"%d",&num_seg);
			fscanf(fpt1,"%s",tempo);
			for(k=0;k<num_seg;k++){
				fscanf(fpt,"%d",&pos);
				fscanf(fpt,"%s",str);
				d=strlen(str);
				for (j=0;j<d;j++){
					ct++;
				}
			}
			if((num_seg!=1)||(ct!=0)){
				i2++;
				snp_ct[i2]=ct+1;
				snp_mat[i2]=malloc(ct*sizeof(int*));
				snp_list[i2]=malloc(ct*sizeof(int*));
				ct=-1;
				for(k=0;k<num_seg;k++){
					fscanf(fpt1,"%d",&pos);
					fscanf(fpt1,"%s",str);
					d=strlen(str);
					for (j=0;j<d;j++){
						ct++;
						snp_list[i2][ct]=pos+j-1;
						col_ct[pos+j-1]++;
						col_ct2++;
						if(str[j]=='1'){
							snp_mat[i2][ct]=1;
						}
						if(str[j]=='2'){
							snp_mat[i2][ct]=2;
						}
						if(str[j]=='3'){
							snp_mat[i2][ct]=3;
						}
						if(str[j]=='4'){
							snp_mat[i2][ct]=4;
						}
					}
				}
				read_ac[i]=0;
				no_acc++;
			}
			else {
				for(k=0;k<num_seg;k++){
					fscanf(fpt1,"%d",&pos);
					fscanf(fpt1,"%s",str);
				}
				read_ac[i]=1;
			}
			fscanf(fpt,"%s",qlty);
			fscanf(fpt1,"%s",qlty);
		}
	}
	fclose(fpt);
	fclose(fpt1);
	int **haplotype;
	haplotype=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		haplotype[iter]=malloc(num_col*sizeof(int));
	}
	for(iter2=0;iter2<num_col;iter2++){
		for (iter=0;iter<K;iter++){
			haplotype[iter][iter2]=1;
		}
	}
	fpt=fopen(argv[1],"r");
	i=-1;i2=-1;
	fscanf(fpt,"%d",&num_read);
	fscanf(fpt,"%d",&num_col);
	int **snp_col;
	snp_col=malloc(num_col*sizeof(int*));
	for(k=0;k<num_col;k++){
		snp_col[k]=malloc(col_ct[k]*sizeof(int));
	}
	for(k=0;k<num_col;k++){
		col_ct[k]=0;
	}
	if(fpt==NULL)
	printf("\nError\n");
	else {
		while(i<num_read-1){
			ct=-1;
			i++;
			fscanf(fpt,"%d",&num_seg);
			fscanf(fpt,"%s",tempo);
			if(read_ac[i]==0){
				i2++;
				for(k=0;k<num_seg;k++){
					fscanf(fpt,"%d",&pos);
					fscanf(fpt,"%s",str);
					d=strlen(str);
					for (j=0;j<d;j++){
						ct++;
						snp_col[pos+j-1][col_ct[pos+j-1]]=i2;
						col_ct[pos+j-1]++;
					}
				}
			}
			else{
				for(k=0;k<num_seg;k++){
					fscanf(fpt,"%d",&pos);
					fscanf(fpt,"%s",str);
					d=strlen(str);
				}
			}
			fscanf(fpt,"%s",qlty);
		}
	}
	fclose(fpt);
	printf("Time to read data %f\n",(double)(clock()-t1)/CLOCKS_PER_SEC);
	num_read=no_acc;
	int *que_read;
	int *que_list;
	int *que_col;
	int iter_cnt;
	int que_ct;
	int *queue;
	int tu, iter3;

	t1=clock();
	row=malloc(10000000*sizeof(int));
	col=malloc(10000000*sizeof(int));
	val=malloc(10000000*sizeof(int));
	queue=malloc(10000000*sizeof(int));

	que_read=malloc(num_read*sizeof(int));
	que_list=malloc(num_read*sizeof(int));
	que_read=malloc(num_read*sizeof(int));
	que_col=malloc(num_col*sizeof(int));
	for(iter=0;iter<num_read;iter++){
		que_read[iter]=0;
	}
	for(iter=0;iter<num_read;iter++){
		que_list[iter]=0;
	}
	iter_cnt=1;
	que_ct=-1;

	int rd_cnt;
	int rd_cnt2;
	int *col_que;
	int colct;
	int *col_arr;
	int *col_arr2;
	int basect;
	int *reord;

	col_que=malloc(num_col*sizeof(int));
	reord=malloc(num_col*sizeof(int));
	col_arr=malloc(num_col*sizeof(int));
	col_arr2=malloc(num_col*sizeof(int));
	for(iter=0;iter<num_col;iter++){
		col_que[iter]=0;
	}
	rd_cnt2=0;
	// File *file_cc_pt; //connected components in dictionary format
	file_cc_pt=fopen("connected_dic.txt","w");

	for(iter=0;iter<num_read;iter++){
		rd_cnt=0;
		basect=0;
		if(que_read[iter]==0){
			que_read[iter]=iter_cnt;
			que_list[iter]=iter_cnt;
			rd_array[rd_cnt]=iter;
			rd_array2[iter]=rd_cnt;
			rd_cnt++;
			basect+=snp_ct[iter];
			for(iter2=0;iter2<snp_ct[iter];iter2++){
				for(iter3=0;iter3<col_ct[snp_list[iter][iter2]];iter3++){
					if(que_list[snp_col[snp_list[iter][iter2]][iter3]]==0){
						que_list[snp_col[snp_list[iter][iter2]][iter3]]=iter_cnt;
						que_ct++;
						queue[que_ct]=snp_col[snp_list[iter][iter2]][iter3];
					}
				}
			}
			while(que_ct!=-1){
				tu=queue[que_ct];
				que_read[tu]=iter_cnt;
				rd_array[rd_cnt]=tu;
				rd_array2[tu]=rd_cnt;
				rd_cnt++;
				basect+=snp_ct[tu];
				que_list[tu]=iter_cnt;
				que_ct--;
				for(iter2=0;iter2<snp_ct[tu];iter2++){
					for(iter3=0;iter3<col_ct[snp_list[tu][iter2]];iter3++){
						if(que_list[snp_col[snp_list[tu][iter2]][iter3]]==0){
							que_list[snp_col[snp_list[tu][iter2]][iter3]]=iter_cnt;
							que_ct++;
							queue[que_ct]=snp_col[snp_list[tu][iter2]][iter3];
						}
					}
				}
			}
			colct=0;
			for(iter2=0;iter2<rd_cnt;iter2++){
				for(iter3=0;iter3<snp_ct[rd_array[iter2]];iter3++){
					if(col_que[snp_list[rd_array[iter2]][iter3]]==0){
						col_que[snp_list[rd_array[iter2]][iter3]]=iter_cnt;
						colct++;
						col_arr[colct]=snp_list[rd_array[iter2]][iter3];
						col_arr2[snp_list[rd_array[iter2]][iter3]]=colct;
						reord[colct-1]=snp_list[rd_array[iter2]][iter3];
					}
				}
			}
			rd_cnt2+=rd_cnt;
			if(colct>1){
				////fprintf("\nNew non-empty block as python dictionary {snp_index: allele}  \n");
				i=0;
				for(iter2=0;iter2<rd_cnt;iter2++){	// iter2-th read in this block
					for(iter3=0;iter3<snp_ct[rd_array[iter2]];iter3++){
						row[i]=iter2+1;	// row[last_i-1]= number of reads
						col[i]=col_arr2[snp_list[rd_array[iter2]][iter3]]; // the snp(column) index  , in special order
                        
						//val[i]=snp_mat[rd_array[iter2]][iter3];  // alleles in the fragment file as an array for this block correspond to the snp index  in col
						i++;
					}
				}
				fprintf(file_cc_pt,"block\n");
				for(int i_frag=0;i_frag<rd_cnt;i_frag++){
					fprintf(file_cc_pt,"%d\n",rd_array[i_frag]+1);
				}


				iter_cnt++;

			} // if(colct>1)
			//else {
			//}
            
            
            
            
            
            
		}// if
	}
	fclose(file_cc_pt); // closing the output of connected componenets

}


