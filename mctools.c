#include <stdlib.h>
#include <stdio.h>
#include "mctools.h"
#include <math.h> 
#include "mymath.h"
#include "ewald_line.h"
#include "io.h"

#define BINS 300
#define ANGLEBINS 100 
#define PI  3.14159265358979323846264338328

int MCstep(double** xyz,double* Qs, int Nions, int Ncat, double stepsize,int Npol, double* box, double** polcoord, double* rions, double* rpol, int*** ioncoord, double* energy)
{

//perform a monte carlo step
double delta_u=0;
int acc=0;
int move_ndx=0;
double x=0,y=0,z=0; //store trial move coordinates here
int ion_type=0;
int temp=0;
//choose move
double move=0
if (Ncat!=0 && Ncat!=Nions)
	move=rand()/(double)RAND_MAX;
if (move<0.5)
{
	transl(stepsize,xyz,&move_ndx, Nions,box,&x,&y,&z);

	if (Qs[move_ndx]>0)
		ion_type=0;
	else
		ion_type=1;

	//check for move overlapping with polymer
	
	if (Npol==1)
		temp=overlap_pol(polcoord[0][0], polcoord[0][1],x,y, rpol[0], rions[ion_type],box);
	else if (Npol==2)
	{
		temp=overlap_pol(polcoord[0][0], polcoord[0][1],x,y, rpol[0], rions[ion_type],box);
		temp=temp+overlap_pol(polcoord[1][0], polcoord[1][1],x,y, rpol[1], rions[ion_type],box);
	}

	temp=temp+overlap_ion(Nions,xyz,move_ndx, rions,Qs,x,y,z,box);

	if (temp>=1)
	{
		
		*energy=0;
		return 0;		
		
	}


	delta_u=calc_delta_e(xyz, Qs,x,y,z, move_ndx);
	acc=accept(delta_u);
	if (acc==1)
	{
		update(xyz,move_ndx,x,y,z,Npol,polcoord,ioncoord, box);
		*energy=delta_u;
	
	}
	else
		*energy=0;
	return acc;
}
else
{
move_ndx=(int)(Ncat*((double)rand() / ((double)(RAND_MAX)+(double)(1)) ));
int move_ndx2=(int)((Nions-Ncat)*((double)rand() / ((double)(RAND_MAX)+(double)(1)) ))+Ncat;
delta_u=calc_deltaE_switch(xyz, Qs, move_ndx,move_ndx2);
//printf(" ndx1 %d q1 %lf  ndx2 %d q2 %lf deltaU %lf \n", move_ndx, Qs[move_ndx],move_ndx2, Qs[move_ndx2], delta_u);
printf("switch");	
	acc=accept(delta_u);
	if (acc==1)
	{
		x=xyz[move_ndx2][0];
		y=xyz[move_ndx2][1];
		z=xyz[move_ndx2][2];
		double x2=xyz[move_ndx][0];
		double y2=xyz[move_ndx][1];
		double z2=xyz[move_ndx][2];
		update(xyz,move_ndx,x,y,z,Npol,polcoord,ioncoord, box);
		update(xyz,move_ndx2,x2,y2,z2,Npol,polcoord,ioncoord, box);
		*energy=delta_u;
	
	}
	else
		*energy=0;
	return acc;
}	


}
//printf("move_ndx %d newkoord %lf %lf %lf delta_e %lf",move_ndx, x,y,z,delta_u);


//-----------------------------------------------------------------
int accept(double delta_u)
{
int acc=0;
double tmp=rand()/(double)RAND_MAX;
//printf("delta u %lf, exp(-deltau) %lf, random %lf\n", delta_u,exp(-delta_u), tmp);
if (exp(-delta_u)>tmp)
	acc=1;
return acc;
}

//-----------------------------------------------------------------
int transl(double step, double** xyz, int* move_ndx, int Nions, double* box, double* x, double* y, double* z)
{
int temp=0;

//translate one aprticle

//choose particle

*move_ndx=(int)(Nions*((double)rand() / ((double)(RAND_MAX)+(double)(1)) ));
//printf("Nions %d attemp move on %d \n",Nions,*move_ndx);
//move
double r[3]={rand()/(double)RAND_MAX, rand()/(double)RAND_MAX, rand()/(double)RAND_MAX};
*x=xyz[*move_ndx][0]+(r[0]-0.5)*step*2;
*y=xyz[*move_ndx][1]+(r[1]-0.5)*step*2;
*z=xyz[*move_ndx][2]+(r[2]-0.5)*step*2;
//*d=5;

put_to_box3D(x, y, z, box);
return temp;

}
//-----------------------------------------------------------------
/*
void pair_transl(double step, double** xyz, int move_ndx, int Nions, int pbc, double* box, double x, double y, double z)
{
//translate a pair of particles



}
*/
//---------------------------------------------------------------

int put_to_box3D(double* x, double* y, double* z, double* box)
{

if (*x<0)
	*x=*x+box[0];	
else if (*x>box[0])
	*x=*x-box[0];

if (*y<0)
	*y=*y+box[1];	
else if (*y>box[1])
	*y=*y-box[1];


if (*z<0)
	*z=*z+box[2];	
else if (*z>box[2])
	*z=*z-box[2];


return 1;
}

//-----------------------------------------------------------------

void update(double** xyz, int move_ndx, double x, double y, double z, int Npol,double** polcoord, int*** ioncoord, double* box)
{
//for now, max distance for ion distribution is the larges box length
//double dx, dy,dz;
//double r2=0;
//double r=0;
//double maxdistance=box[1];
//double angle=0;
/*
for (int i=1;i<3;i++)
{
	if (maxdistance<box[i])
		maxdistance=box[i];

}

for (int i=0;i<Npol;i++)
{
//Note, distance not corrected with hard sphere radius, convenient for ploting
// Angle wrt positive x-axis (1,0,0)
	dx = x - polcoord[i][0]; dx -= dround(dx/box[0])*box[0];
	dy = y - polcoord[i][1]; dy -= dround(dy/box[1])*box[1]; 
	dz = 0;
	r2 = SQR(dx) + SQR(dy) + SQR(dz);
	r=sqrt(r2);
	//going through quadrats of coordinate system, scaling acos 0,pi->0,2pi
	angle=acos(dx/r);

	if (dx>0 && dy<0)
		angle=2*PI-angle;
	else if (dx<0 && dy<0)
		angle=(PI-angle)+PI;

	
	ioncoord[i][move_ndx][0]=(int)(angle/(2*PI)*ANGLEBINS);
	ioncoord[i][move_ndx][1]=(int)(r*BINS/maxdistance);

} */
xyz[move_ndx][0]=x;
xyz[move_ndx][1]=y;
xyz[move_ndx][2]=z;




}

//-----------------------------------------------------------------
int update_iondist(int Nions, int Npol, double* Qs, int*** ioncoord, int*** cationdist, int*** aniondist,int*** cationdistsum, int*** aniondistsum, int*** aniondistsum2, int*** cationdistsum2, int*** totiondistsum, int*** totiondistsum2,  int set2zero)
{
if (set2zero==1)
{
  for (int i=0;i<Npol; i++)
  {  	
		for (int k=0;k<ANGLEBINS;k++)
		{
			for (int l=0;l<BINS;l++)
			{
				
					cationdist[i][k][l]=0;
					aniondist[i][k][l]=0;
					cationdistsum[i][k][l]=0;
					aniondistsum[i][k][l]=0;
					cationdistsum2[i][k][l]=0;
					aniondistsum2[i][k][l]=0;
					totiondistsum[i][k][l]=0;
					totiondistsum2[i][k][l]=0;

			}
		}
			
  }

return 1;

}	
else{
for (int i=0;i<Npol; i++)
{  
	for (int j=0;j<Nions;j++)
	{
	if (Qs[j]>0)
		cationdist[i][ioncoord[i][j][0]][ioncoord[i][j][1]]++;
	else
		aniondist[i][ioncoord[i][j][0]][ioncoord[i][j][1]]++;
		
	}	
}
  for (int i=0;i<Npol; i++)
  {  	
		for (int k=0;k<ANGLEBINS;k++)
		{
			for (int l=0;l<BINS;l++)
			{
				cationdistsum[i][k][l]=cationdistsum[i][k][l]+cationdist[i][k][l];
				cationdistsum2[i][k][l]=cationdistsum2[i][k][l]+cationdist[i][k][l]*cationdist[i][k][l];
				aniondistsum[i][k][l]=aniondistsum[i][k][l]+aniondist[i][k][l];
				aniondistsum2[i][k][l]=aniondistsum2[i][k][l]+aniondist[i][k][l]*aniondist[i][k][l];
				totiondistsum[i][k][l]=totiondistsum[i][k][l]+cationdist[i][k][l]-aniondist[i][k][l];
				totiondistsum2[i][k][l]=totiondistsum2[i][k][l]+(cationdist[i][k][l]-aniondist[i][k][l])*(cationdist[i][k][l]-aniondist[i][k][l]);
				cationdist[i][k][l]=0;
				aniondist[i][k][l]=0;
			}
		}
			
  }
return 0;

}


}
//--------------------------------------------------------------------
int update_ionpress(int Npol, int Nions,double* Qs, double** xyz, double* rpol, double** polcoord, double* rions,int*** press_extrap,double* box, int setzero)
{
double dr[]={0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08, 0.085, 0.09, 0.095, 0.1,0.105,0.11,0.115,0.12, 0.125,0.130};
double dx=0;
double dy=0;
double r2=0;
double r=0;
double angle=0;
if (setzero==1)
{

for (int i=0;i<Npol;i++)
	for (int j=0;j<ANGLEBINS;j++)
		for (int k=0;k<50;k++)	
			press_extrap[i][j][k]=0;

return 1;
}
else
{
	for (int i=0;i<Npol;i++)
	{
	//Note, distance not corrected with hard sphere radius, convenient for ploting
	// Angle wrt positive x-axis (1,0,0)
		for (int j=0;j<Nions; j++)
		{
		dx = xyz[j][0] - polcoord[i][0]; dx -= dround(dx/box[0])*box[0];
		dy = xyz[j][1] - polcoord[i][1]; dy -= dround(dy/box[1])*box[1]; 
		
		r2 = SQR(dx) + SQR(dy);
		r=sqrt(r2);
		//going through quadrats of coordinate system, scaling acos 0,pi->0,2pi
		angle=acos(dx/r);

		if (dx>0 && dy<0)
			angle=2*PI-angle;
		else if (dx<0 && dy<0)
			angle=(PI-angle)+PI;
	
			for (int k=0;k<25;k++)
			{
			if (Qs[j]>0 && r<rions[0]+rpol[i]+dr[k])
				press_extrap[i][(int)(angle/(2*PI)*ANGLEBINS)][k]++;
			if (Qs[j]<0 && r<rions[1]+rpol[i]+dr[k])
				press_extrap[i][(int)(angle/(2*PI)*ANGLEBINS)][k+25]++;
			}		
		}


	}

return 0;
}

}



//--------------------------------------------------------------------
double calc_delta_e(double** xyz, double* Qs, double x, double y, double z, int move_ndx)
{
// calc ewald energy+other neccesary energies

double new_coord[3]={x,y,z};
double old_coord[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
      double Er1=Ewald_r_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx); 
      double Ek1=Ewald_k_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx);
//      double Ed1=Ewald_dipol(xyz, Qs, 0,0,0);
//     double Old_e=calc_energy(xyz,Qs);
//double Old_er=Ewald_r_space_line(xyz, Qs,0,0,0);
//double Old_ek=Ewald_k_space_line(xyz, Qs,0,0,0);		
xyz[move_ndx][0]=new_coord[0];
xyz[move_ndx][1]=new_coord[1];
xyz[move_ndx][2]=new_coord[2];
      double Er2=Ewald_r_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
      double Ek2=Ewald_k_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);      
//      double Ed2=Ewald_dipol(xyz, Qs, 0,0,0);
//    double New_e=calc_energy(xyz,Qs);
//double New_er=Ewald_r_space_line(xyz, Qs,0,0,0);
//double New_ek=Ewald_k_space_line(xyz, Qs,0,0,0);	
//printf("old_e %lf new_e %lf delta_e %lf partial_new %lf part_old %lf part_delta %lf\n", Old_e, New_e, New_e-Old_e,(Er2+Ek2),(Er1+Ek1), (Er2+Ek2)-(Er1+Ek1) );   
//printf("1 real %lf k %lf real_tot %lf k_tot %lf\n", Er1, Ek1, Old_er, Old_ek);
//printf("2 real %lf k %lf real_tot %lf k_tot %lf\n", Er2, Ek2, New_er, New_ek);
xyz[move_ndx][0]=old_coord[0];
xyz[move_ndx][1]=old_coord[1];
xyz[move_ndx][2]=old_coord[2];


//return (Er2+Ek2+Ed2)-(Er1+Ek1+Ed1);
//dipole off
return (Er2+Ek2)-(Er1+Ek1);
}
//-----------------------------------------------------
double calc_deltaE_switch(double** xyz, double* Qs, int move_ndx, int move_ndx2 )
{
// calc ewald energy+other neccesary energies
double old_coord[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
double old_coord2[3]={xyz[move_ndx2][0],xyz[move_ndx2][1],xyz[move_ndx2][2]};
double new_coord[3]={xyz[move_ndx2][0],xyz[move_ndx2][1],xyz[move_ndx2][2]};
double new_coord2[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
//energy in old conf for particle1
//      double Er11=Ewald_r_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx); 
//      double Ek11=Ewald_k_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx);
//energy in old conf for particle2
//      double Er12=Ewald_r_space_part_line(xyz, Qs, old_coord2, Qs[move_ndx2], move_ndx2); 
//      double Ek12=Ewald_k_space_part_line(xyz, Qs, old_coord2, Qs[move_ndx2], move_ndx2);
//      double Ed1=Ewald_dipol(xyz, Qs, 0,0,0);
//     double Old_e=calc_energy(xyz,Qs);	
double dm_k1=Ewald_k_space_doublemove(xyz, Qs, old_coord, Qs[move_ndx],old_coord2, Qs[move_ndx2]);
double dm_r1=Ewald_r_space_doublemove(xyz, Qs, old_coord, Qs[move_ndx], move_ndx,old_coord2,Qs[move_ndx2], move_ndx2);
	
xyz[move_ndx][0]=new_coord[0];
xyz[move_ndx][1]=new_coord[1];
xyz[move_ndx][2]=new_coord[2];
xyz[move_ndx2][0]=new_coord2[0];
xyz[move_ndx2][1]=new_coord2[1];
xyz[move_ndx2][2]=new_coord2[2];

//      double Er21=Ewald_r_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
//      double Ek21=Ewald_k_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
//     double Er22=Ewald_r_space_part_line(xyz, Qs, new_coord2, Qs[move_ndx2], move_ndx2);
//      double Ek22=Ewald_k_space_part_line(xyz, Qs, new_coord2, Qs[move_ndx2], move_ndx2);        
//      double Ed2=Ewald_dipol(xyz, Qs, 0,0,0);
//    double New_e=calc_energy(xyz,Qs);
double dm_k2=Ewald_k_space_doublemove(xyz, Qs, new_coord, Qs[move_ndx],new_coord2, Qs[move_ndx2]);
double dm_r2=Ewald_r_space_doublemove(xyz, Qs, new_coord, Qs[move_ndx], move_ndx,new_coord2,Qs[move_ndx2], move_ndx2);	
//printf("old_e %lf new_e %lf delta_e %lf partial_new %lf part_old %lf part_delta %lf\n", Old_e, New_e, New_e-Old_e,(Er2+Ek2),(Er1+Ek1), (Er2+Ek2)-(Er1+Ek1) );   
//printf("real %lf k %lf dipole %lf\n", Er2-Er1, Ek2- Ek1, Ed2-Ed1);
xyz[move_ndx][0]=old_coord[0];
xyz[move_ndx][1]=old_coord[1];
xyz[move_ndx][2]=old_coord[2];
xyz[move_ndx2][0]=old_coord2[0];
xyz[move_ndx2][1]=old_coord2[1];
xyz[move_ndx2][2]=old_coord2[2];
//printf("%lf %lf %lf %lf\n", dm_r1, dm_r2,dm_k1, dm_k2);

//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",(dm_r2+dm_k2)-(dm_k1+dm_r1),(Er21+Ek21+Er22+Ek22)-(Er11+Ek11+Er12+Ek12),New_e-Old_e);
//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",dm_r2-dm_r1,(Er21+Er22)-(Er11+Er12),New_e-Old_e);
//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",dm_k2-dm_k1,(Ek21+Ek22)-(Ek11+Ek12),New_e-Old_e);
//return (Er2+Ek2+Ed2)-(Er1+Ek1+Ed1);
//dipole off
return (dm_r2+dm_k2)-(dm_k1+dm_r1);
}
//--------------------------------------------------------------

double calc_energy(double** xyz, double* Qs)
{
// calc ewald energy+other neccesary energies
double er=0;
double ek=0;
double ed=0;
er=Ewald_r_space_line(xyz, Qs,0,0,0);
ek=Ewald_k_space_line(xyz, Qs,0,0,0);
//ed=Ewald_dipol(xyz,Qs,0,0,0);
//dipole off;
ed=0;
//printf("total energy er %lf ek %lf ed %lf\n",er,ek,ed);
return er+ek+ed;
}
//--------------------------------------------------------------
void force_on_polymer(double** xyz, double* Qs, double  force[][3] , int Npol)
{

double Fxpi_r=0;
double Fxpi_k=0;

double Fypi_r=0;
double Fypi_k=0;

double Fxpp_r=0;
double Fxpp_k=0;

double Fypp_r=0;
double Fypp_k=0;

for (int l=0;l<Npol;l++)
{

	ion_to_cylinder_r(xyz, Qs, l,&Fxpi_r, &Fypi_r);
	ion_to_cylinder_k(xyz, Qs,l,&Fxpi_k, &Fypi_k);


	force[l][0]=Fxpi_r+Fxpi_k;	
	force[l][1]=Fypi_r+Fypi_k;
//printf("pol %d Fxr %lf Fxk, %lf Fyr %lf Fyk %lf \n", l,Fxpi_r, Fxpi_k, Fypi_r, Fypi_k);
//printf("pol %d F_x %lf F_y %lf \n", l,force[l][0], force[l][1]);
 Fxpi_r=0;
 Fxpi_k=0;

 Fypi_r=0;
 Fypi_k=0;	

	if (Npol==2)
	{	
	 Fxpp_r=0;
	 Fxpp_k=0;

	 Fypp_r=0;
	 Fypp_k=0;
		cylinder_to_cylinder_k(l, &Fxpp_k, &Fypp_k);
		cylinder_to_cylinder_r(&Fxpp_r, &Fypp_r,l);
		force[l][0]=force[l][0]+Fxpp_r+Fxpp_k;	
		force[l][1]=force[l][1]+Fypp_r+Fypp_k;
	
//	printf("pol %d Fxppr %lf Fxkpp, %lf Fyrpp %lf Fykpp %lf \n",l, Fxpp_r, Fxpp_k, Fypp_r, Fypp_k);
//	printf("pol %d F_x %lf F_y %lf \n", l,force[l][0], force[l][1]);
	}
}



}

//---------------------------------------------------------------

int overlap_pol(double px, double py, double x, double y, double rpol, double r_ion, double* box)
{
double dx=x-px;
double dy=y-py;

dx -= dround(dx/box[0])*box[0];
dx -= dround(dx/box[1])*box[1];

int temp=0;

if (sqrt((dy)*(dy)+(dx)*(dx))<(r_ion+rpol))
	temp=1;

return temp;
}

//------------------------------------------------------------------

int overlap_ion(int Np, double** xyz, int move_ndx, double* rions, double* Qs, double x, double y, double z, double* box)
{

int moved_ion_type;
int ion_type;

double dx;
double dy;
double dz;

/*printf("box %lf %lf %lf\n", box[0],box[1],box[2]);
printf("rions %lf %lf\n", rions[0],rions[1]);
printf("Np %d", Np);*/

if (Qs[move_ndx]>0)
	moved_ion_type=0;
else
	moved_ion_type=1;

for (int i=0;i<Np;i++)
{
		
	if (Qs[i]>0)
		ion_type=0;
	else
		ion_type=1;
	//printf("distance %lf limit %lf\n", SQR(xyz[i][0]-xyz[move_ndx][0])+SQR(xyz[i][1]-xyz[move_ndx][1])+SQR(xyz[i][2]-xyz[move_ndx][2]),SQR(rions[ion_type]+rions[moved_ion_type]));
	//printf("x %lf y %lf z %lf\n", xyz[i][0], xyz[i][1],xyz[i][2]);
	dx=xyz[i][0]-x;
	dy=xyz[i][1]-y;
	dz=xyz[i][2]-z;

	dx -= dround(dx/box[0])*box[0];
	dy -= dround(dy/box[1])*box[1];
	dz -= dround(dz/box[2])*box[2];

	if (move_ndx!=i && (SQR(dx)+SQR(dy)+SQR(dz)<SQR(rions[ion_type]+rions[moved_ion_type])))
	{
		return 1;
	} 	
		
}
return 0;
}
//----------------------------------------------------------------------------

void gen_polymerconf(double dist_pol, int Npol, int Ncat, int Nan, double** polcoord, double* box, double* Qs, double* q_iontype, double* rion, double* rpol, double* taus, double** xyz)
{
double x=0;
double y=0;
double z=0;
int q=0;
int Ncat_t=Ncat;
int Nan_t=Nan;	
int k=0;
int Nions=Ncat+Nan;

	if (Npol==2)
	{


		polcoord[0][0]=box[0]/2.0+dist_pol/2.0;
		polcoord[0][1]=box[1]/2;

		polcoord[1][0]=box[0]/2.0-dist_pol/2.0;
		polcoord[1][1]=box[1]/2;		

	}
	if (Npol==1)
	{
		polcoord[0][0]=box[0]/2;
		polcoord[0][1]=box[1]/2;
	
	}
	

	
	
	while (k<Nions)
	{
		// randomise for coordinate directon
	
		x=(rand()/(double)RAND_MAX*box[0]);
		y=(rand()/(double)RAND_MAX*box[1]);
		z=(rand()/(double)RAND_MAX*box[2]);
		
		// randomise for cation/anion 0/1
		
		if (Ncat_t>0)
		{
			q=q_iontype[0];
			
		}				
		else
		{
			q=q_iontype[1];
			
		}
		
		
		//check for overlaps, talle on helpompikin toteutus	
		if (q==q_iontype[1])
		{	
			Qs[k]=q_iontype[1];
			xyz[k][0]=x;
			xyz[k][1]=y;
			xyz[k][2]=z;
			if (overlap_pol(polcoord[0][0],polcoord[0][1],x,y,rpol[0],rion[1],box)==1)
				continue;		
			if (Npol>1)
			{
			 if (overlap_pol(polcoord[1][0], polcoord[1][1],x,y,rpol[1],rion[1],box)==1) 
				continue;
			}
			if(overlap_ion(Nions,xyz,k,rion, Qs,x,y,z,box))
				continue;
				
			k++;
			Nan_t--;		
		}
		else
		{
			Qs[k]=q_iontype[0];
			xyz[k][0]=x;
			xyz[k][1]=y;
			xyz[k][2]=z;
			

			if (overlap_pol(polcoord[0][0],polcoord[0][1],x,y,rpol[0],rion[0],box)==1)
				continue;	
			if (Npol>1)
			{
				if( overlap_pol(polcoord[1][0], polcoord[1][1],x,y,rpol[1],rion[0], box)==1)
				continue;
			}
			if(overlap_ion(Nions,xyz,k,rion, Qs,x,y,z,box))
				continue;			
			
			k++;
			Ncat_t--;
			
		}
						

	} 	
printf("random configuration generated\n");
printf("%d polymers\n", Npol);
printf("%d cations and %d anions\n", Ncat-Ncat_t, Nan-Nan_t);
printf("box %lf %lf %lf taus %lf %lf q_iontype %lf %lf \n",box[0], box[1], box[2], taus[0], taus[1], q_iontype[0], q_iontype[1]); 
if (Npol>1)
	printf("total charge %f\n", box[2]*taus[0]+box[2]*taus[1]+Ncat*q_iontype[0]+Nan*q_iontype[1]);
else
	printf("total charge %f\n", box[2]*taus[0]+Ncat*q_iontype[0]+Nan*q_iontype[1]);

write_initial_gro("generated_configuration.gro",xyz, Nions, Npol, polcoord, box, Qs);



}


