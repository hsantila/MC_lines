#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include <errno.h>


#define BINS 300
#define ANGLEBINS 100
#define PI  3.14159265358979323846264338328

//------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------------------------

int read_input(char* filein,char* fileout,int* relax,int* nsteps, int* step_upd, int* gro_outp, int* var_outp, double* stepsize, double* alpha, double* rcut, int* kcut, int* Npol, int* Ncat, int* Nan, double* rion, double* qion, double* rpol, double* dist_pol, double* taus, double* box)
{


FILE *inp;
FILE *outp;
int ret=0;
char dataid[500];
double dataval=0.0;

taus[0]=0;
taus[1]=0;
qion[0]=0;
qion[1]=0;

if ((inp = fopen(filein,"r")) != NULL)
{
	
	while(ret!=EOF)
	{
		ret=fscanf(inp,"%s", dataid);
		//ei tunnista # merkkiä kommentiksi, huono implementaatio.
		if (strstr(dataid,"Npolymers")!=NULL)
		{
			
			ret=fscanf(inp, " %lf\n", &dataval);
			*Npol=(int)dataval;
			
		}
		if(strstr(dataid,"Nanion")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*Nan=(int)dataval;
		}
		if(strstr(dataid,"Ncation")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*Ncat=(int)dataval;
		}
		if(strstr(dataid,"qcation")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			qion[0]=dataval;
		}
		if(strstr(dataid,"qanion")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			qion[1]=dataval;
		}
		if(strstr(dataid,"r_anion")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			rion[1]=dataval;
		}
		if(strstr(dataid,"r_cation")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			rion[0]=dataval;
		}
		if(strstr(dataid,"r_pol1")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			rpol[0]=dataval;
		}
		if(strstr(dataid,"r_pol2")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			rpol[1]=dataval;
		}
		if(strstr(dataid,"taupol1")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			taus[0]=dataval;
		}
		if(strstr(dataid,"taupol2")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			taus[1]=dataval;
		}
		if(strstr(dataid,"dist_poly")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*dist_pol=dataval;
		}
		if(strstr(dataid,"BOXX")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[0]=dataval;
		}
		if(strstr(dataid,"BOXY")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[1]=dataval;
		}
		if(strstr(dataid,"BOXZ")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[2]=dataval;
		}
		if(strstr(dataid,"stepsize")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*stepsize=dataval;
		}
		if(strstr(dataid,"stepupdateinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*step_upd=(int)dataval;
		}
		if(strstr(dataid,"dataprintoutinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*var_outp=(int)dataval;
		}
		if(strstr(dataid,"groprintoutinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*gro_outp=(int)dataval;
		}
		if(strstr(dataid,"RELAXATIONSTEPS")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*relax=(int)dataval;
		}
		if(strstr(dataid,"MCsteps")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*nsteps=(int)dataval;
		}
		if(strstr(dataid,"Ealpha")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			printf("alpha %lf\n", dataval);
			*alpha=dataval;
		}
		if(strstr(dataid,"Rcut")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*rcut=dataval;
		}
		if(strstr(dataid,"EwaldKx")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[0]=(int)dataval;
		}
		if(strstr(dataid,"EwaldKy")!=NULL)
		{
			
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[1]=(int)dataval;
		}
		if(strstr(dataid,"EwaldKz")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[2]=(int)dataval;
		}
//		if(strstr(dataid,""))
//			=(int)dataval;

	}




	fclose(inp);
	
outp=fopen(fileout,"w");
 fprintf(outp,"# BOX_X\nBOXX %lf\n\n",  box[0]);
 fprintf(outp,"# BOX_Y\nBOXY %lf\n\n",  box[1]);
 fprintf(outp,"# BOX_Z\nBOXZ %lf\n\n",  box[2]);
 fprintf(outp,"# NPolymer\nNpolymers %d\n\n",  *Npol);
 fprintf(outp,"# line charge 1\ntaupol1 %lf\n\n",  taus[0]);
 fprintf(outp,"# line charge 2\ntaupol2 %lf\n\n",taus[1] );
 fprintf(outp,"# Anion charge\nqanion %lf\n\n",qion[1]);
 fprintf(outp,"# Cation charge\nqcation %lf\n\n",qion[0]);
 fprintf(outp,"# Amount of cations\nNcation %d\n\n",*Ncat);
 fprintf(outp,"# Amount of anios\nNanion %d\n\n", *Nan);
 fprintf(outp,"# Anion radius\nr_cation %lf\n\n",rion[1]);
 fprintf(outp,"# Cation radius\nr_anion %lf\n\n",rion[0]);
 fprintf(outp,"# Polymer 1 radius\nr_pol1 %lf\n\n",rpol[0]);
 fprintf(outp,"# Polymer 2 radius\nr_pol2 %lf\n\n",rpol[1]);
 fprintf(outp,"# Polymer distance\ndist_poly %lf\n\n",*dist_pol);
 fprintf(outp,"# MC steps\nMCsteps %d\n\n", *nsteps );
 fprintf(outp,"# GRO printout interval \ngroprintoutinterval %d\n\n", *gro_outp);
 fprintf(outp,"# Data printout interval \ndataprintoutinterval %d\n\n", *var_outp);
 fprintf(outp,"# Step update interval \nstepupdateinterval %d\n\n", *step_upd);
 fprintf(outp,"# Initial steps to be disregarded \nRELAXATIONSTEPS %d\n\n", *relax);
 fprintf(outp,"# stepsize\nstepsize %lf\n\n", *stepsize );
 fprintf(outp,"# K vectos in x\nEwaldKx %d\n\n",kcut[0]);
 fprintf(outp,"# K vectos in y\nEwaldKy %d\n\n",kcut[1]);
 fprintf(outp,"# K vectos in z\nEwaldKz %d\n\n",kcut[2]);
 fprintf(outp,"# Direct space cutoff\nRcut %lf\n\n",*rcut);
 fprintf(outp,"# Ewald summation alpha\nEalpha %lf\n\n",*alpha );
// fprintf(outp,"# \n %d\n\n", );
fclose(outp);

	return 1;
}
else
{
	printf("Cannot read input parameters\n");
	return 0;
}

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------

int write_initial_gro(char* file, double** xyz, int Nions, int Npol, double** polcoord, double* box, double* Qs)
{
FILE *gro;

if ((gro = fopen(file,"w")) != NULL)
{


int res=100;
int atomnumber=1;
int resnumber=1;
char* str;
char atomn[5];
char atomtype[5];
double zkoord=0;


fprintf(gro,"Polymer(s) and ions, initial configuration\n");
fprintf(gro, "%d\n", res*Npol+Nions);

//write polymers
for (int i=0;i<Npol;i++)
{
	strcpy(atomtype, "P");
	sprintf(atomn,"%d",i+1);
	strcat(atomtype, atomn);
	for (int j=1+res*i;j<=res+i*res;j++)
	{	
		
		fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "POL",atomtype, j, polcoord[i][0], polcoord[i][1], zkoord+box[2]/res*(j-1-res*i));
		atomnumber=atomnumber+1;

	}
	
}
resnumber=Npol+1;

//write ions

for (int i=0;i<Nions;i++)
{
	if (Qs[i]>0)
		str="CAT";
	else
		str="ANI";
	fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resnumber, str,str,atomnumber,xyz[i][0], xyz[i][1],xyz[i][2]);	
	resnumber=resnumber+1;
	atomnumber=atomnumber+1;
}

//write box

fprintf(gro,"%f %f %f\n", box[0], box[1], box[2]);
fclose(gro);
return 1;
}
else
{
	printf("Cannot open the gro file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

int write_gro(char* file,int step, double** xyz, int Nions, int Npol, double** polcoord, double* box, double* Qs)
{
FILE *gro;

if ((gro = fopen(file,"a")) != NULL)
{


int res=100;
int atomnumber=1;
int resnumber=0;
char* str;
char atomn[5];
char atomtype[5];
double zkoord=0;


fprintf(gro,"Polymer(s) and ions, t= %d\n", step);
fprintf(gro, "%d\n", res*Npol+Nions);

//write polymers
for (int i=0;i<Npol;i++)
{
	strcpy(atomtype, "P");
	sprintf(atomn,"%d",i+1);
	strcat(atomtype, atomn);
	for (int j=1+res*i;j<=res+i*res;j++)
	{	
		
		fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "POL",atomtype, j, polcoord[i][0], polcoord[i][1], zkoord+box[2]/res*(j-1-res*i));
		
		atomnumber=atomnumber+1;

	}
	
}
resnumber=Npol+1;

//write ions

for (int i=0;i<Nions;i++)
{
	if (Qs[i]>0)
		str="CAT";
	else
		str="ANI";
	fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resnumber, str,str,atomnumber,xyz[i][0], xyz[i][1],xyz[i][2]);	
	resnumber=resnumber+1;
	atomnumber=atomnumber+1;
}

//write box

fprintf(gro,"%f %f %f\n", box[0], box[1], box[2]);
fclose(gro);
return 1;
}
else
{
	printf("Cannot open the gro file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
int read_gro(char* file, double** xyz, double** polcoord, double* q_iontype, double* Qs, int* Npol, int* Ncat, int* Nan, double* box)
{
int maxlen=100;
char line[maxlen];
FILE *gro;
int Natoms=0;
int resnumber=0;
char resname[5];
int atomnumber=0;
char atomname[5];
int tmp=0;
int k=0;
if ((gro = fopen(file,"r")) != NULL)
{
fgets(line, maxlen, gro);
fscanf(gro,"%d\n",&Natoms);
	
	for (int i=0;i<Natoms;i++)
	{
		fscanf(gro, "%5d%5s %5s %5d %lf %lf %lf\n",&resnumber,resname, atomname,&atomnumber,&xyz[k][0],&xyz[k][1],&xyz[k][2]);
		
		if (strstr(resname, "POL")!=NULL)
		{
			if (resnumber==1)
			{
				polcoord[0][0]=xyz[k][0];
				polcoord[0][1]=xyz[k][1];
				tmp=1;
			}
			else
			{
				polcoord[1][0]=xyz[k][0];
				polcoord[1][1]=xyz[k][1];
				tmp=2;
			}
			
			continue;
		}
		if (strstr(atomname, "ANI")!=NULL)
		{
			Qs[k]=q_iontype[1];	
			*Nan=*Nan+1;
			k++;
		}
		if (strstr(atomname, "CAT")!=NULL)
		{
			Qs[k]=q_iontype[0];
			*Ncat=*Ncat+1;
			k++;
		}
	}
*Npol=tmp;
fscanf(gro, "%lf %lf %lf", &box[0],&box[1],&box[2]);
printf("Read %d polymers, %d cations and %d anions\n", *Npol,*Ncat,*Nan);
fclose(gro);
return 1;
}
else
{
	printf("Cannot read from gro file, The reason *may* be %s\n", strerror(errno));
	return 0;
}
}

/*--------------------------------------------------------------------------------------------------------------------*/


int write_energy(char* file,int step, double* energy)
{
FILE *ene;
if ((ene = fopen(file,"a")) != NULL)
{
	fprintf(ene,"%d %lf\n", step, *energy);
fclose(ene);
return 1;
}
else
	printf("Cannot open the energy file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}

/*-------------------------------------------------------------------------------------------------------------------*/

int write_stepsize(char* file,int step, double stepsize)
{
FILE *ene;
if ((ene = fopen(file,"a")) != NULL)
{
	fprintf(ene,"%d %lf\n", step, stepsize);
fclose(ene);
return 1;
}
else
	printf("Cannot open the stepsize file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}

/*--------------------------------------------------------------------------------------------------------------------*/

int write_acc(char* file,int step, double acc)
{
FILE *ene;
if ((ene = fopen(file,"a")) != NULL)
{
	fprintf(ene,"%d %lf\n", step, acc);
fclose(ene);
return 1;
}
else
	printf("Cannot open the acceptance ration file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------*/


int write_iondist(char* file_base,int Nave, int Npol, int*** iondistsum, double* box)
{
FILE *dist;
double maxdistance=box[1];
char full_filename[60];
char tmpstr[20];

for (int i=1;i<3;i++)
{
	if (maxdistance<box[i])
		maxdistance=box[i];

}
for (int p=0;p<Npol;p++)
{
	sprintf(tmpstr, "pol%d_", p);	
	strcpy(full_filename, tmpstr);
	strcat(full_filename, file_base);
	if ((dist = fopen(full_filename,"a")) != NULL)
	{
		fprintf(dist,"%d", Nave);

		for (int l=0;l<BINS;l++)
			fprintf(dist,"\t%lf",(l+0.5)/BINS*(double)maxdistance);

		fprintf(dist,"\n ");

			for (int j=0;j<ANGLEBINS;j++)
			{
				fprintf(dist,"%lf", 2.0*PI*(j+0.5)/(double)ANGLEBINS);
						
				for (int k=0;k<BINS;k++)
				{
				
					fprintf(dist,"\t%lf", (double)iondistsum[p][j][k]/((double)Nave));
				}
				fprintf(dist, "\n");
			}

	fprintf(dist, "\n");		
	fclose(dist);
	}
	else
	{	printf("Cannot open the ion distribution file for writing, The reason *may* be %s\n", strerror(errno));
		return 0;
	}

}
return 1;
}
//-----------------------------------------------------------------------------------------------------------------------
int write_pressextrp(int Npol,char* file_base, int*** press_extra, int avrg_count)
{
FILE *press;

double angle=0;

char tmpstr[20];
char full_filename[60];
for (int i=0;i<Npol;i++)
{
sprintf(tmpstr, "pol%d_", i);	
strcpy(full_filename, tmpstr);
strcat(full_filename, file_base);

	if ((press= fopen(full_filename,"a")) != NULL)
	{
		for (int j=0;j<ANGLEBINS;j++)
		{
			fprintf(press,"%lf\t", angle=2.0*PI*(j+0.5)/(double)ANGLEBINS);
			for (int k=0;k<50;k++)					
				fprintf(press,"  %lf",(double)press_extra[i][j][k]/(double)avrg_count);
		fprintf(press,"\n");	
		}
	fclose(press);
	}
	else
	{
	printf("Cannot open the acceptance ration file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
	}	

}	


return 1;
}

/*------------------------------------------------------------------------------------------------------------------------*/
int write_forces(char* file, int step, double forces[][3])
{
static int firstwrite=0;
	
FILE *ene;

if ((ene = fopen(file,"a")) != NULL)
{
	if (firstwrite==0)
		fprintf(ene,"# step \t Fx_pol1 \t Fy_pol1 \t Fz_pol1 \t Ftot_pol1 \t Fx_pol2 \t Fy_pol2 \t Fz_pol2 \t Ftot_pol2\n");
fprintf(ene,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", step, forces[0][0], forces[0][1],forces[0][2],sqrt(forces[0][0]*forces[0][0]+forces[0][1]*forces[0][1]+forces[0][2]*forces[0][2]),forces[1][0], forces[1][1],forces[1][2],sqrt(forces[1][0]*forces[1][0]+forces[1][1]*forces[1][1]+forces[1][2]*forces[1][2]));

fclose(ene);
firstwrite+=1;


return 1;
}
else
	printf("Cannot open the acceptance ration file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;

}
