#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "string.h"
#include "io.h"
#include "mctools.h"
#include "ewald_line.h"
#include <time.h>

#define BINS 300
#define ANGLEBINS 100 // 60
#define PI  3.14159265358979323846264338328

int main(int argc, char *argv[])
{
//step parameters
int nstep=0;	//number of steps
double stepsize=0; 
int start_t=0;
int relax=0;
int acc=0;
int acc_sum=0;
//int lb;

//output&update intervals

int step_upd=0;
int gro_outp=0;
int var_outp=0;

//filenames

char* file_in;
char* file_out_base;
char* file_gro;
char file_out_gro[40];
char file_out_energ[40];
char file_out_ionAn[40];
char file_out_ionCat[40];
char file_out_totions[40];
char file_out_ionAn2[40];
char file_out_ionCat2[40];
char file_out_totions2[40];
char file_out_forces[40];
char file_out_extr[40];

//indexes and indicators for input

int is_g=0;
int is_i=0;
int is_o=0;
int is_e=0;
int is_t=0;


int ndx_g=0;
int ndx_i=0;
int ndx_o=0;
int ndx_e=0;

//int ndx_t=0;



//Ewald parameters
double alpha=0;
//double coulomb=0;
double rcut=0;
int kcut[3];
//isotropic ewald convergence
int kmax_c=0;
int kmin_c=0;

//system spesification


int Nions=0;
int Nan=0;
int Ncat=0;
double q_iontype[2];
double rion[2];

int Npol=0;
double taus[2];
double rpol[2];
double** polcoord;
double dist_pol;

double box[3];
double** xyz;
double* Qs; 
//periodicity on/off
//int pbc=0;



//other random tool variables
//double rnd=0;
int temp=0;
int seed=1;

//analysis variables
double energy=0;
double etot=0;
int*** ioncoord;
int*** aniondist;
int*** cationdist;
int*** aniondistsum;
int*** cationdistsum;
int*** aniondistsum2;
int*** cationdistsum2;
int*** totiondistsum;
int*** totiondistsum2;
int*** press_extra;

int avg_count=0;
double forces[2][3]={{0,0,0},{0,0,0}};
//------------------------------------------------------------------------

//                Reading input parameters

//------------------------------------------------------------------------



for (int i=1; i<argc; i++){
    if (strcmp(argv[i],"-g")==0){
      is_g = 1;
      ndx_g = i;
    }
    if (strcmp(argv[i],"-i")==0){
      is_i = 1;
      ndx_i = i;
    }
    if (strcmp(argv[i],"-o")==0){
      is_o = 1;
      ndx_o = i;
    }
    if (strcmp(argv[i],"-e")==0){
	  is_e = 1;
     	 ndx_e = i;
    }
    if (strcmp(argv[i],"-s")==0){
	nstep =(int)strtod(argv[i+1], NULL); // If parameter file contains steps, this is overwritten
    }
    if (strcmp(argv[i],"-r")==0){
	seed =(int)strtod(argv[i+1], NULL); // optional seed for random number generator
    }		
    if (strcmp(argv[i],"-t")==0){
	is_t=1;
	start_t =(int)strtod(argv[i+1], NULL); // If gro file given
    }
    	
}

//check for input options
if (is_o==0 || is_i==0){
	printf("Invalid input\n");
	printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e  [Rcut] [kmin] [kmax]\n", argv[0] );
	
	exit(0);
}

if((is_g==0 && is_t==1) || (is_g==1 && is_t==0))
{
	printf("Invalid input in optional parameters\n");
	printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e [Rcut] [kmin] [kmax]\n", argv[0] );
	
	exit(0);
}


if(is_i==1)
{
	file_in=argv[ndx_i+1];
}
if(is_o==1)
{
	file_out_base=argv[ndx_o+1];
}
if(is_g==1)
{
	file_gro=argv[ndx_g+1];
}

strcpy(file_out_gro, file_out_base);
strcat(file_out_gro,".gro");
strcpy(file_out_energ, file_out_base);
strcat(file_out_energ,".enrg");
strcpy(file_out_ionAn, file_out_base);
strcat(file_out_ionAn,".An");
strcpy(file_out_ionCat, file_out_base);
strcat(file_out_ionCat,".Cat");
strcpy(file_out_ionAn2, file_out_base);
strcat(file_out_ionAn2,".An2");
strcpy(file_out_ionCat2, file_out_base);
strcat(file_out_ionCat2,".Cat2");
strcpy(file_out_totions, file_out_base);
strcat(file_out_totions,".totion");
strcpy(file_out_totions2, file_out_base);
strcat(file_out_totions2,".totion2");
strcpy(file_out_forces,file_out_base);
strcat(file_out_forces,".force");
strcpy(file_out_extr,file_out_base);
strcat(file_out_extr,".extrp");
//------------------------------------------------------------------------------------------

//               Reading information from files, memory allocation for dynamic and system dependent matrixes

//------------------------------------------------------------------------------------------



//for the moment, max 2 polymers
polcoord=malloc(sizeof(double*)*2);
	for (int i=0; i<2; i++)
		polcoord[i]=malloc(sizeof(double)*2);
polcoord[0][0]=0;
polcoord[1][0]=0;
polcoord[0][1]=0;
polcoord[1][1]=0;		
	

temp=read_input(file_in,file_out_base,&relax,&nstep, &step_upd, &gro_outp, &var_outp, &stepsize, &alpha, &rcut, kcut, &Npol, &Ncat, &Nan, rion, q_iontype, rpol,&dist_pol,taus, box);

Nions=Ncat+Nan;
if(temp==0)
{
	printf("Cannot read input file\n");
	exit(0);
}

//tee failsafe muuttujille

cationdist=malloc(sizeof(int**)*Npol);
aniondist=malloc(sizeof(int**)*Npol);
cationdistsum=malloc(sizeof(int**)*Npol);
aniondistsum=malloc(sizeof(int**)*Npol);
cationdistsum2=malloc(sizeof(int**)*Npol);
aniondistsum2=malloc(sizeof(int**)*Npol);
totiondistsum=malloc(sizeof(int**)*Npol);
totiondistsum2=malloc(sizeof(int**)*Npol);
press_extra=malloc(sizeof(int**)*Npol);
	for (int j=0;j<Npol;j++)
	{
		cationdist[j]=malloc(sizeof(int*)*ANGLEBINS);
		aniondist[j]=malloc(sizeof(int*)*ANGLEBINS);
		cationdistsum[j]=malloc(sizeof(int*)*ANGLEBINS);
		aniondistsum[j]=malloc(sizeof(int*)*ANGLEBINS);
		cationdistsum2[j]=malloc(sizeof(int*)*ANGLEBINS);
		aniondistsum2[j]=malloc(sizeof(int*)*ANGLEBINS);
		totiondistsum[j]=malloc(sizeof(int*)*ANGLEBINS);
		totiondistsum2[j]=malloc(sizeof(int*)*ANGLEBINS);
		press_extra[j]=malloc(sizeof(int*)*ANGLEBINS);
		for (int i=0;i<ANGLEBINS;i++)
		{
			cationdist[j][i]=malloc(sizeof(int)*BINS);
			aniondist[j][i]=malloc(sizeof(int)*BINS);
			cationdistsum[j][i]=malloc(sizeof(int)*BINS);
			aniondistsum[j][i]=malloc(sizeof(int)*BINS);
			cationdistsum2[j][i]=malloc(sizeof(int)*BINS);
			aniondistsum2[j][i]=malloc(sizeof(int)*BINS);
			totiondistsum2[j][i]=malloc(sizeof(int)*BINS);
			totiondistsum[j][i]=malloc(sizeof(int)*BINS);
			press_extra[j][i]=malloc(sizeof(int)*50);

		}
	}	
xyz=malloc(sizeof(double*)*Nions);
	for (int i=0; i<Nions;i++)
		xyz[i]=calloc(1,sizeof(double)*3);

Qs=calloc(1,sizeof(double)*Nions);


ioncoord=malloc(sizeof(int**)*Npol);
	for (int j=0;j<Npol;j++)
	{
		ioncoord[j]=malloc(sizeof(int*)*Nions);	
			for (int i=0;i<Nions;i++)
				ioncoord[j][i]=malloc(sizeof(int)*3);
	}

if (is_g==1)
{
//readin configuration
int Npol_gro=0;
int Ncat_gro=0;
int Nan_gro=0;
read_gro(file_gro, xyz, polcoord, q_iontype, Qs, &Npol_gro, &Ncat_gro, &Nan_gro, box);

	if(Npol!=Npol_gro || Ncat!=Ncat_gro || Nan!=Nan_gro)
	{

	printf("Missmatch in amount of polymers or ions in input file and geometry read from gro\n");
	exit(0);
	}
printf("read in configuration\n");
if (Npol==0)
printf("system total charge %lf\n", Ncat*q_iontype[0]+Nan*q_iontype[1]);
if (Npol==1)
printf("system total charge %lf\n", Ncat*q_iontype[0]+Nan*q_iontype[1]+taus[0]*box[2]);
if (Npol==2)
printf("system total charge %lf\n", Ncat*q_iontype[0]+Nan*q_iontype[1]+(taus[0]+taus[1])*box[2]);
for (int i=0;i<Nions;i++)
{
	printf("%d %lf %lf %lf %lf\n", i, xyz[i][0],xyz[i][1],xyz[i][2], Qs[i]);	

}

}
else
{
//generate configuration
gen_polymerconf(dist_pol, Npol, Ncat, Nan, polcoord, box, Qs, q_iontype,rion, rpol, taus, xyz);

}
//------------------------------------------------------------------------

//               Ewald convergence test

//------------------------------------------------------------------------

if(is_e==1)
{
char file_abs[50];
char file_move[50];
char file_abs_parts[50];
char tmpstr[25];





/*FILE* move=fopen("move_energy.txt","w");
FILE* abs_ene=fopen("total_energy","w");
fprintf(move,"#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9\n");
fprintf(abs_ene,"#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9\n");*/
FILE* abs;
FILE* move;
FILE* abs_parts;


double e1=0;
double delta_e=0;
double delta_e1=0;
double delta_e2=0;
double delta_e3=0;




double ek_ii=0;
double ek_pi=0;
double ek_pp=0;
double ek_self=0;
double er_ii=0;
double er_pi=0;
double er_pp=0;
double ek_tot=0;
double er_tot=0;

int move_ndx;

double amin=0;
double amax=0;
double alpha=0;
char direct;
int k_const=0;
int karr[3];
//int ndx=0;
move_ndx=Nions-1;



	if (ndx_e+8<argc)
	{
		printf("Not enough input parameters for ewald convergence test\n");
		printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e [Rcut] [amin] [amax] [kmin] [kmax] [direction x, y or z] [kconst]\n", argv[0] );
		exit(0);
		//kconst=0 symmetry
	}	 	
	else
	{

		rcut = strtod(argv[ndx_e+1], NULL);
		amin=(double)strtod(argv[ndx_e+2], NULL);
		amax=(double)strtod(argv[ndx_e+3], NULL);
		kmin_c = (int)strtod(argv[ndx_e+4], NULL);
		kmax_c = (int)strtod(argv[ndx_e+5], NULL);
		direct=argv[ndx_e+6][0];

		if (direct=='x')
			ndx=0;
		if (direct=='y')
			ndx=1;
		if (direct=='z')
			ndx=2;

		k_const = (int)strtod(argv[ndx_e+7], NULL);
		printf("rcut %lf amin %lf amax %lf kmin %d kmax %d\n", rcut, amin, amax, kmin_c, kmax_c);	
		//perform ewald convergence test		
		for (int i=kmin_c; i<kmax_c+1; i++)
		{
			//fprintf(move,"%d", i);
			//fprintf(abs_ene,"%d", i);

			strcpy(file_abs,"tot_energy_conv");
			strcpy(file_move,"move_energy_conv");
			strcpy(file_abs_parts,"tot_energy_parts");

			sprintf(tmpstr,"_Rcut%1.2lfk%d%c%d",rcut,k_const,direct,i);
		
			strcat(file_move,tmpstr);
			strcat(file_move,".xvg");



			
			strcat(file_abs,tmpstr);
			strcat(file_abs,".xvg");
			strcat(file_abs_parts,tmpstr);
			strcat(file_abs_parts,".xvg");		

			printf("filuabs %s\n", file_abs);
			printf("filumove %s\n", file_move);


			abs=fopen(file_abs,"w");
			move=fopen(file_move,"w");
			abs_parts=fopen(file_abs_parts,"w");

			for (int j=0;j<amax*100;j++)
			{
				alpha = amin+(j+1)/100.0;
				if (k_const==0)
				{
				karr[0]=i;
				karr[1]=i;
				karr[2]=i;
				}
				else
				{

				karr[0]=i;
				karr[1]=i;
				karr[2]=k_const;
			
				}
				
	
				printf("kx %d ky %d kz %d\n", karr[0],karr[1], karr[2]);

				Ewald_init_line(box, Nions,Npol, alpha,rcut,karr , 1, 1, taus, polcoord);

				e1=calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]+0.1;
				xyz[move_ndx][1]=xyz[move_ndx][1]+0.5;
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.5;
	
				delta_e=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-0.1;
				xyz[move_ndx][1]=xyz[move_ndx][1]-0.5;
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.5;


				xyz[move_ndx][0]=xyz[move_ndx][0]+0.5;
				xyz[move_ndx][1]=xyz[move_ndx][1]+0.5;
	
	
				delta_e1=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-0.5;
				xyz[move_ndx][1]=xyz[move_ndx][1]-0.5;
	

	
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.5;
	
				delta_e2=e1-calc_energy(xyz, Qs);
			
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.5;

				xyz[move_ndx][0]=xyz[move_ndx][0]+7;
				xyz[move_ndx][1]=xyz[move_ndx][1]+7;
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.1;
	
				delta_e3=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-7;
				xyz[move_ndx][1]=xyz[move_ndx][1]-7;
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.1;

				fprintf(move,"%.12lf %.12lf %.12lf %.12lf %.12lf\n", alpha, delta_e, delta_e1, delta_e2, delta_e3);	
				fprintf(abs,"%lf %lf\n",alpha, e1);

				//ek_tot=Ewald_r_space(xyz, Qs, 0);
				ek_tot=Ewald_k_space_line_parts(xyz, Qs, &ek_ii, &ek_pi,  &ek_pp, &ek_self);
				er_tot=Ewald_r_space_line_parts(xyz, Qs, &er_ii, &er_pi, &er_pp);
				fprintf(abs_parts,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", alpha, e1,ek_tot, er_tot,ek_pi, ek_ii, ek_pp, er_pi, er_ii, er_pp, ek_self);
			
	
										
			}			

				//fprintf(move,"\n");
				//fprintf(abs_ene,"\n");
			fclose(move);
			fclose(abs);
			fclose(abs_parts);		
		}
		exit(0);	
	}


}

//------------------------------------------------------------------------

//              Do monte carlo!

//------------------------------------------------------------------------

//initiate ioncoords
for (int i=0;i<Nions;i++)
{
if (Qs[i]>0)
	temp=0;
else
	temp=1;
	update(xyz,i, xyz[i][0], xyz[i][1], xyz[i][2], Npol,polcoord, ioncoord, box);
}


//set iondistribution to zero
//update_iondist(Nions, Npol, Qs, ioncoord, cationdist, aniondist,cationdistsum, aniondistsum, aniondistsum2, cationdistsum2, totiondistsum, totiondistsum2, 1);
update_ionpress(Npol, Nions, Qs, xyz, rpol, polcoord, rion,press_extra,box, 1);	
//initiate ewald
Ewald_init_line(box, Nions,Npol, alpha,rcut, kcut, 1, 1, taus, polcoord);

//initioate total energy
etot=calc_energy(xyz,Qs);

printf("askel %d energia %lf\n",-1, etot);
forces[0][0]=0;
forces[0][1]=0;
forces[0][2]=0;
forces[1][0]=0;	
forces[1][1]=0;
forces[1][2]=0;	
force_on_polymer(xyz, Qs, forces, Npol);
temp=write_forces(file_out_forces, -1, forces);	
forces[0][0]=0;
forces[0][1]=0;
forces[0][2]=0;
forces[1][0]=0;	
forces[1][1]=0;
forces[1][2]=0;	

//set the random number generator to runtime, doesn't work with supercomputers
//srand ( (unsigned)time ( NULL ) );
//set the random number generator to input seed
srand(seed);

//THE monte carlo loop
for (int i=start_t; i<=nstep;i++)
{

	acc=MCstep(xyz,Qs,Nions,Ncat, stepsize,Npol, box, polcoord, rion, rpol,ioncoord, &energy);
	acc_sum=acc_sum+acc;
	etot=etot+energy;
		
	
	//if (acc==1)
//	{
		
		//printf("askel %d etot %lf deltaE %lf systemtot %lf \n",i, etot,energy,calc_energy(xyz,Qs));
			
//	}
	//else
	//	printf(" move rejected \n");
	//printf("step %d acc_sum %d\n", i, acc_sum);
	
	if(i<relax)
	{
		//do relaxation
		
		continue;

	}

	

	if((i+1)%step_upd == 0)
	{
		//update stepsize
		if((1.0*acc_sum/step_upd)>0.75)
		{
			stepsize=stepsize*10.0;
			if (stepsize>box[0] && stepsize>box[1] && stepsize>box[2])
			{
				printf("Attempt to increase step size over box size on timestep %d\n",i);
				stepsize=stepsize/10.0;
			}
			else 
				printf("Stepsize updated to %f on timestep %d\n", stepsize,i);			
							
		}
		if((1.0*acc_sum/step_upd)<0.25)
		{
						
			stepsize=stepsize/10.0;
			if(stepsize<0.1)
			{
				printf("Attempt to decrease steps size under limit 0.1 on timestep %d\n",i);
				stepsize=stepsize*10.0;
				
			}
			else
				printf("Stepsize updated to %f on timestep %d\n", stepsize,i);	

		}
		acc_sum=0.0;
	}
	
	if (i%gro_outp==0)
	{
	//basic printouts, energy, gro and force	
		temp=write_gro(file_out_gro,i, xyz, Nions, Npol, polcoord, box, Qs);
		if (temp==0)
			exit(0);
		printf("step %d acc_sum %d energy %f\n", i, acc_sum, etot);
		temp=write_energy(file_out_energ,i, &etot);
		if (temp==0)
			exit(0);
		force_on_polymer(xyz, Qs, forces, Npol);
		temp=write_forces(file_out_forces, i, forces);
		forces[0][0]=0;
		forces[0][1]=0;
		forces[0][2]=0;
		forces[1][0]=0;
		forces[1][1]=0;
		forces[1][2]=0;	
		if (temp==0)
			exit(EXIT_FAILURE);
		
	}
	if ((i)%100==0)
	{
	 //collecting averages
	 //  hard-coded to 100, if var_outp aka dataprintoutinterval<100 produces zeors and nans	 	
		/*update_iondist(Nions, Npol, Qs, ioncoord, cationdist, aniondist,cationdistsum, aniondistsum, aniondistsum2, cationdistsum2, totiondistsum, totiondistsum2, 0);*/
		update_ionpress(Npol, Nions, Qs, xyz, rpol, polcoord, rion,press_extra,box, 0);	
		avg_count=avg_count+1;
	}		
	if(i%var_outp==0)
	{
	//add analysis tool printouts
	//for now, just iondist, keep var_output sparse, othervise huge files 
		
		temp=write_pressextrp(Npol,file_out_extr, press_extra, avg_count);
		if (temp==0)
			exit(EXIT_FAILURE);

	/*	temp=write_iondist(file_out_ionAn, avg_count,Npol, aniondistsum, box);
		temp=write_iondist(file_out_ionAn2, avg_count,Npol, aniondistsum2, box);
		temp=write_iondist(file_out_ionCat,avg_count, Npol, cationdistsum, box);
		temp=write_iondist(file_out_ionCat2,avg_count, Npol, cationdistsum2, box);
		temp=write_iondist(file_out_totions, avg_count,Npol, totiondistsum, box);
		temp=write_iondist(file_out_totions2, avg_count,Npol, totiondistsum2, box);
		//iondist to zero
		update_iondist(Nions, Npol, Qs, ioncoord, cationdist, aniondist,cationdistsum, aniondistsum, aniondistsum2, cationdistsum2, totiondistsum, totiondistsum2, 1);*/
		update_ionpress(Npol, Nions, Qs, xyz, rpol, polcoord, rion,press_extra,box, 1);	
		avg_count=0;
		if (temp==0)
			exit(EXIT_FAILURE);
		


			
	}

		

}
//--------------------------------------------------------------------------------

//				Memory handling

//--------------------------------------------------------------------------------


	for (int j=0;j<Npol;j++)
	{

		for (int i=0;i<ANGLEBINS;i++)
		{
			free(cationdist[j][i]);
			free(aniondist[j][i]);
			free(cationdistsum[j][i]);
			free(aniondistsum[j][i]);
			free(cationdistsum2[j][i]);
			free(aniondistsum2[j][i]);
			free(totiondistsum[j][i]);
			free(totiondistsum2[j][i]);
			free(press_extra[j][i]);

		}
		
		free(cationdist[j]);
		free(aniondist[j]);
		free(cationdistsum[j]);
		free(aniondistsum[j]);
		free(cationdistsum2[j]);
		free(aniondistsum2[j]);
		free(totiondistsum[j]);
		free(totiondistsum2[j]);
		free(press_extra[j]);

	}

free(cationdist);
free(aniondist);
free(cationdistsum);
free(aniondistsum);
free(cationdistsum2);
free(aniondistsum2);
free(totiondistsum);
free(totiondistsum2);
free(press_extra);

for (int i=0; i<2;i++)
	free(polcoord[i]);
free(polcoord);


for (int i=0; i<Nions;i++)
	free(xyz[i]);
free(xyz);

for (int j=0;j<Npol;j++)
{
	for (int i=0; i<Nions;i++)
		free(ioncoord[j][i]);
	free(ioncoord[j]);
}
free(ioncoord);

free(Qs);




}
