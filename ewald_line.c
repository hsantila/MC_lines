/*------------------------------------------------------------

  SUBUNIT:  ewald.c

  NEEDS:    ewald.h, mymath.h, mymath.c

  PURPOSE:  Puts a simple calculation of the real- and reciprocal
            space contribution to the Ewald sum to the disposal
	    of the user.

  COMMENTS: - Not time-optimized!
            - Contains a net charge correction in case the 
	      system is not electrostatically neutral.
	    - if forces should be computed, the three addresses of 
	      the force arrays have to be given to the routines
	      Ewald_r_space, Ewald_k_space and  Ewald_dipol. 
	      If one does NOT want this, one MUST instead give
	      three 0 to the routines - this tells them to
	      skip over the force computation.

 
  CONTENT:  Ewald_init :
              Initialize the Ewald sum, store some internal constants
	      and precompute the influence function.

	    Ewald_r_space : 
              Compute the real space contribution 

            Ewald_k_space :
              Compute the reciprocal space contribution

	    Ewald_dipol :
	      Compute the dipol correction

  USE:      The four functions mentioned above have the following calling syntax :
  
            void Ewald_init(double [length of system], 
                            int    [number of particles], 
			    double [Ewald parameter alpha],
			    double [real-space-cutoff], 
			    int    [reciprocal-space-cutoff], 
			    double [prefactor of the Coulomb potential, e.g.
			            bjerrumlength*temperature],
			    double [epsilon of surrounding dielectric])

	    double Ewald_r_space(double *[array of x-coordinates of the particles], 
                                 double *[array of y-coordinates of the particles], 
				 double *[array of z-coordinates of the particles], 
				 double *[array of particle charges], 
				 double *[array of x-component of the forces], 
				 double *[array of y-component of the forces], 
				 double *[array of y-component of the forces]);
		   The function returns the real space 
		   contribution to the Coulomb energy.

	    double Ewald_k_space(double *[array of x-coordinates of the particles], 
                                 double *[array of y-coordinates of the particles], 
				 double *[array of z-coordinates of the particles], 
				 double *[array of particle charges], 
				 double *[array of x-component of the forces], 
				 double *[array of y-component of the forces], 
				 double *[array of y-component of the forces]);
		   The function returns the reciprocal space 
		   contribution to the Coulomb energy, including
		   self energy and net charge correction.

	    double Ewald_k_subset(double *[array of x-coordinates of the particles], 
	                          double *[array of y-coordinates of the particles], 
				  double *[array of z-coordinates of the particles], 
				  double *[array of particle charges], 
				  int [start-number for subset], 
				  int [end-number for subset]);
		   The function calculates the contribution to the
		   Ewald k-space energy (including self energy and
		   net charge correction) resulting from interactions
		   of charges with index in {start,start+1,...,end}.
		   Note that 'start' and 'end' are INCLUDED in the list!

	    double Ewald_dipol(double *[array of x-coordinates of the particles], 
	                       double *[array of y-coordinates of the particles], 
			       double *[array of z-coordinates of the particles], 
			       double *[array of particle charges], 
			       double *[array of x-component of the forces], 
			       double *[array of y-component of the forces], 
			       double *[array of y-component of the forces]);
		   The function returns the dipol 
		   correction to the Coulomb energy.

------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ewald_line.h"
#include "mymath.h"

/*----------------------------------------------------------------------*/

/* Some variables, internal to ewald.c, which are shared among the routines: */

static double  L3;    /* box lengths, its inverse and volume */
static int     NP;         /* number of particles */
static int     Npol;	   /*number of cylinders*/
static double  alpha;      /* Ewald parameter */
static double  rmax,rmax2, kmax2; /* real space cutoff and its square */
static double  prefactor;  /* prefactor of the Coulomb potential */
static double  sum_q_2;    /* SQUARE of SUM of charges */
static double  sum_q2;     /* SUM of SQUARE of charges */
static double  dipfac;     /* 2*PI / ((1 + 2*epsilon) * L^3) */

/* stores the influence function: */
static double  Ghat[Maxkmax+1][Maxkmax+1][Maxkmax+1];

static int kmax[3]; //unisotropic cutoff matrixes
static double Ls[3];
 
/* needed for the Fourier transforms: */
static double  e_re[2*Maxkmax+1][2*Maxkmax+1][2*Maxkmax+1];
static double  e_im[2*Maxkmax+1][2*Maxkmax+1][2*Maxkmax+1];

/*Cylinder line charge spesification*/
static double pxy[2][2]; 	/*cylinder coordinates*/
static double taus[2];	/*cylinder line charges*/

/*----------------------------------------------------------------------*/

static void calculate_influence_function(void);

/*----------------------------------------------------------------------*/

void Ewald_init_line(double* length, int particlenumber,int polymern, double Ewaldalpha,
		double realspacecutoff, int* kcut, 
		double coulombprefactor, double epsilon, double* line_c, double** ccoord) 
{
  static double PI = 3.14159265358979323846264338328;

  Ls[0]        = length[0];
  Ls[1]		= length[1];
  Ls[2]	=length[2];


  L3        = Ls[0]*Ls[1]*Ls[2];
  NP        = particlenumber;
  Npol	    = polymern;
  alpha     = Ewaldalpha;
  rmax      = realspacecutoff;
  rmax2     = SQR(rmax);
kmax[0]= kcut[0];
kmax[1]=kcut[1];
kmax[2]=kcut[2];
//for pseudo psherical cutoff
  kmax2     = (SQR(kmax[0])+SQR(kmax[1])+SQR(kmax[2]))/3;
  prefactor = coulombprefactor;
  dipfac    = 2.0 * PI / ( (1.0 + 2.0*epsilon) * L3 );
taus[0]=line_c[0];
taus[1]=line_c[1];

pxy[0][0]=ccoord[0][0];
pxy[0][1]=ccoord[0][1];
pxy[1][0]=ccoord[1][0];
pxy[1][1]=ccoord[1][1];


  if (kmax[0] > Maxkmax && kmax[1]>Maxkmax && kmax[2]>Maxkmax) {
    fprintf(stderr,"Subroutine 'Ewald-init' complains:\n"
	    "Desired maximal k-vector (%d %d %d) is larger\n"
	    "than what the program can cope with (%d)\n"
	    "Program terminated!\n\n",kmax[0], kmax[1], kmax[2],Maxkmax);
    exit(1);
  }
    
  calculate_influence_function();
  
  fprintf(stderr,"'Ewald' successfully initialized!\n");
	
}
/* internal function */

static void calculate_influence_function(void);

/*----------------------------------------------------------------------*/

/* external functions */

void   Ewald_init(double length, int particlenumber, double Ewaldalpha,
		  double realspacecutoff, int reciprocalspacecutoff, 
		  double coulombprefactor, double epsilon);

double Ewald_r_space_line(double** xyz, double *Q,
		     double *Fx, double *Fy, double *Fz);

double Ewald_k_space_line(double** xyz, double *Q,
		     double *Fx, double *Fy, double *Fz);

double Ewald_dipol(double** xyz, double *Q,
		   double *Fx, double *Fy, double *Fz);

/*----------------------------------------------------------------------*/

static void calculate_influence_function(void)
{
  /* calculates the influence function 2.0/L^2 * exp(-(PI*n/(alpha*L))^2)/n^2
     as a function of lattice vector n (NOT k=2*PI*n/L). This is stored in the
     array Ghat[Maxkmax+1][Maxkmax+1][Maxkmax+1]. For symmetry reasons only
     one octant is actually needed. */
  
  int    nx,ny,nz,ii=0;
  double n_2;

  static double PI = 3.14159265358979323846264338328;
  
  fprintf(stderr,"calculating influence function for the standard Ewald method using\n"
	         "L=%lf %lf %lf\tkmax=%d %d %d\talpha=%lf ...",Ls[0] ,Ls[1], Ls[2],kmax[0], kmax[1], kmax[2],alpha);

  for (nx = 0; nx <= kmax[0]; nx++)
    for (ny = 0; ny <= kmax[1]; ny++)
      for (nz = 0; nz <= kmax[2]; nz++) {
	if ((nx==0) && (ny==0) && (nz==0)) Ghat[nx][ny][nz] = 0.0;
	else {
	  n_2 = SQR(2*PI)*(SQR(nx/Ls[0]) + SQR(ny/Ls[1]) + SQR(nz/Ls[2]));
	  Ghat[nx][ny][nz] = prefactor * exp(-n_2/(4*SQR(alpha))) / n_2;
	
	  ii=ii+1;
	}
      }
  fprintf(stderr,"\n");
  printf("\n wavevectrors%le\n",1.0*ii);
}

/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

double Ewald_dipol(double** xyz, double *Q, 
		   double *Fx, double *Fy, double *Fz)
{
  double dipx=0.0, dipy=0.0, dipz=0.0;
  double d1=prefactor*dipfac, d2=0.0;
  int    i;

  for (i=0; i<NP; i++) {
   // printf("ind %d \n",i);
   // printf("QS[i] %lf\n", Qs[i]);
   // printf("x %lf y %lf  z %lf\n", xyz[i][0], xyz[i][1], xyz[i][2]);			
    dipx += Q[i]*xyz[i][0];
    dipy += Q[i]*xyz[i][1];
    dipz += Q[i]*xyz[i][2];
  }

  d2 += d1 * ( SQR(dipx) + SQR(dipy) + SQR(dipz) );

  if (Fx||Fy||Fz) {   /* if forces should be computed... */
    d1   *= 2.0;
    dipx *= d1;
    dipy *= d1;
    dipz *= d1;
    for (i=0; i<NP; i++) {
      Fx[i] -= Q[i] * dipx;
      Fy[i] -= Q[i] * dipy;
      Fz[i] -= Q[i] * dipz;
    }
  }
  return d2;
}
/*----------------------------------------------------------------------*/

double Ewald_k_space_line(double** xyz, double *Q, 
		     double *Fx, double *Fy, double *Fz)
{
  int    i;
  int    nx,ny,nz;
  double s,h,c,re,im;
  double d1,d2=0.0,d3=0.0,ci=0.0, cc=0.0;
  //double ci2=0;
  double hpol[2];	
  static double PI   = 3.14159265358979323846264338328;
  static double wupi = 1.77245385090551602729816748334;
   
  d1 = 2*PI/L3;

  sum_q_2 = sum_q2 = 0.0;


  for (i=0; i<NP; i++) {
    sum_q_2 += Q[i];
    sum_q2  += SQR(Q[i]);
  }
  sum_q_2 *= sum_q_2;

//  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
//    for (ny=-kmax[1]; ny<=kmax[1]; ny++)
//      for (nz=-kmax[2]; nz<=kmax[2]; nz++)
//	e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = 0.0;

  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++)
      for (nz=-kmax[2]; nz<=kmax[2]; nz++) {
	e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = 0.0;
		hpol[0]=2.0 * PI * (nx*pxy[0][0]/Ls[0] + ny*pxy[0][1]/Ls[1]);
 	 	//if (Npol==2) removed if for vectorization, pxy sould be initialized
		hpol[1]=2.0 * PI * (nx*pxy[1][0]/Ls[0] + ny*pxy[1][1]/Ls[1]);	
	for (i=0; i<NP; i++)
	//remove this if, for nonsymmetrik 
	//  if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	  {
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]);
	    e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] += Q[i]*cos(h);
	    e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] -= Q[i]*sin(h);
	    d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];

	    if (nz==0)
	    {	
		for (int l=0; l<Npol;l++)
		{	
	    	//ci2=ci2+(e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*cos(hpol[l])- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*sin(hpol[l]))*taus[l]*Ls[2]*2*d2;
		ci=ci+d2*taus[l]*Ls[2]*cos(-hpol[l]+h)*2*Q[i];	
	    	}

	  }
	}
	for (i=0; i<NP; i++)
	//remove this if, for nonsymmetrik k
	// if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	{
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]) ;
	    

	    c = cos(h);
	    s = sin(h);

	    d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
	    re = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c 
		- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s) * d2;
	

	    if (Fx||Fy||Fz) {   /* if forces should be computed... */
	      im = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s 
                  + e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c) * d2;

	      d2 = im * Q[i];

	      Fx[i] += 2*PI*nx * d2/Ls[0];
	      Fy[i] += 2*PI*ny * d2/Ls[1];
	      Fz[i] += 2*PI*nz * d2/Ls[2];
	    }
	   

	    d3 += Q[i] * re;
	 // cylinder-cylinder
	}
	if(nz==0)
	{
          
	//if for pseudo symmetric boundary
	// if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2)
	//{
	//	for (int l=0; l<Npol;l++)
	//	{	
	    	//cc=cc+taus[l]*taus[l]*Ls[2]*Ls[2]*d2;
	//	}
		if (Npol==2)
			cc=cc+2*taus[0]*taus[1]*Ls[2]*Ls[2]*cos(hpol[0]-hpol[1])*d2;
	//}
	}	
      }

//printf("k-space d3 %lf cc %lf ci %lf cc %lf\n", d1*d3,d1*cc,d1*ci, d1); 
  d3=d3+cc+ci;

  
  d3 *= d1;

  /* self energy and net charge correction _added_ not substractes as normal. Taus should be assigned zero if no value */
  d3 -= prefactor * ( sum_q2 * alpha / wupi  -  (SQR(taus[1]*Ls[2])+SQR(taus[0]*Ls[2])) * PI / (2.0*L3*SQR(alpha)) );
// no self correction, no net charge

//with no net charge
//  d3-= prefactor * ( sum_q2 * alpha / wupi );
// with no self correction
//d3-=prefactor*(sum_q_2 * PI / (2.0*L3*SQR(alpha)));
// with all
//d3 -= prefactor * ( sum_q2 * alpha / wupi  +  sum_q_2 * PI / (2.0*L3*SQR(alpha)) );

  /*  printf("\n%le", sum_q2); */
  
  return d3;
}

/*----------------------------------------------------------------------*/
double Ewald_k_space_line_parts(double** xyz, double *Q, 
		     double* e_ii, double* e_pi, double* e_pp, double* e_self)
{

//separates different kontributions in k-space
  int    i;
  int    nx,ny,nz;
  double s,h,c,re;
  double d1,d2=0.0,d3=0.0,ci=0.0, cc=0.0;
//double ci2=0;
  double hpol[2];	
  static double PI   = 3.14159265358979323846264338328;
  static double wupi = 1.77245385090551602729816748334;
   
  d1 = 2*PI/L3;

  sum_q_2 = sum_q2 = 0.0;


  for (i=0; i<NP; i++) {
    sum_q_2 += Q[i];
    sum_q2  += SQR(Q[i]);
  }
  sum_q_2 *= sum_q_2;



  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++)
      for (nz=-kmax[2]; nz<=kmax[2]; nz++) {
	e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = 0.0;
		hpol[0]=2.0 * PI * (nx*pxy[0][0]/Ls[0] + ny*pxy[0][1]/Ls[1]);
 	 	//if (Npol==2) removed if for vectorization, pxy sould be initialized
		hpol[1]=2.0 * PI * (nx*pxy[1][0]/Ls[0] + ny*pxy[1][1]/Ls[1]);	
	for (i=0; i<NP; i++)
	//remove this if, for nonsymmetrik 
	 // if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	  {
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]);
	    e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] += Q[i]*cos(h);
	    e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] -= Q[i]*sin(h);
	    d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];

	    if (nz==0)
	    {	
		for (int l=0; l<Npol;l++)
		{	
	    	//ci2=ci2+(e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*cos(hpol[l])- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*sin(hpol[l]))*taus[l]*Ls[2]*2*d2;
		ci=ci+d2*taus[l]*Ls[2]*cos(-hpol[l]+h)*2*Q[i];	
	    }

	  }
	}
	for (i=0; i<NP; i++)
	//remove this if, for nonsymmetrik k
	// if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	{
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]) ;
	    

	    c = cos(h);
	    s = sin(h);

	    d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
	    re = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c 
		- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s) * d2;
	


	
	   

	    d3 += Q[i] * re;
	 // cylinder-cylinder
	}
	if(nz==0)
	{
           	//if for pseudo symmetric boundary
	// if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2)
	//{
	//	for (int l=0; l<Npol;l++)
	//	{	
	    	//cc=cc+taus[l]*taus[l]*Ls[2]*Ls[2]*d2;
	//	}
		if (Npol==2)
			cc=cc+2*taus[0]*taus[1]*Ls[2]*Ls[2]*cos(hpol[0]-hpol[1])*d2;
	//}
	

	}	
      }

//printf("k-space d3 %lf cc %lf ci %lf cc %lf\n", d1*d3,d1*cc,d1*ci, d1*ci2); 

*e_ii=d1*d3;
*e_pi=d1*ci;
*e_pp=d1*cc;
*e_self=-prefactor * ( sum_q2 * alpha / wupi );
  d3=d3+cc+ci;

  
  d3 *= d1;

  /* self energy and net charge correction note the signs and if no value assigned tau should be zero */
  d3 -= prefactor * ( sum_q2 * alpha / wupi  -  (SQR(taus[1]*Ls[2])+SQR(taus[0]*Ls[2]))* PI / (2.0*L3*SQR(alpha)) );
// no self correction, no net charge

//with no net charge
//  d3-= prefactor * ( sum_q2 * alpha / wupi );
// with no self correction
//d3-=prefactor*(sum_q_2 * PI / (2.0*L3*SQR(alpha)));
// with all
//d3 -= prefactor * ( sum_q2 * alpha / wupi  +  sum_q_2 * PI / (2.0*L3*SQR(alpha)) );

  /*  printf("\n%le", sum_q2); */


  return d3;

}
/*----------------------------------------------------------------------*/
double Ewald_r_space_line_parts(double** xyz, double *Q, 
		     double* e_ii, double* e_pi, double* e_pp)
{
  int    t1,t2;
  double dx,dy,dz,r,r2,ar, dpix, dpiy;
  double d1,d2,d3=0.0, d3_ci=0, d3_cc=0;
  double rpi2=0;
  

  for (t1 = 0; t1 < NP-1; t1++)   /* Quick and Dirty N^2 loop */
  {
 //cylinder-ion contribution
    for (int l=0;l<Npol;l++)
    {	
	dpix=pxy[l][0]-xyz[t1][0]; dpix-=dround(dpix/Ls[0])*Ls[0];
	dpiy=pxy[l][1]-xyz[t1][1]; dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));
	//e1(x)=-ei(-x)
	if (rpi2<=rmax2)
   		d3_ci-= Q[t1]*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));
	
	
    }
    	

    for (t2 = t1 + 1; t2<NP; t2++) {
      dx = xyz[t1][0] - xyz[t2][0]; dx -= dround(dx/Ls[0])*Ls[0];
      dy = xyz[t1][1] - xyz[t2][1]; dy -= dround(dy/Ls[1])*Ls[1]; 
      dz = xyz[t1][2] - xyz[t2][2]; dz -= dround(dz/Ls[2])*Ls[2];
      r2 = SQR(dx) + SQR(dy) + SQR(dz);
      /*      printf("%d %d %le \n",t1+1,t2+1,sqrt(r2));*/
      /*      printf("%d  %d  %le  %le\n",t1+1,t2+1,r2,rmax2);*/

	if (r2 <= rmax2) {
	  //printf("Forces");
	  ar = alpha * (r = sqrt(r2));
	  d3 += ( d2 = (d1 = Q[t1] * Q[t2] * prefactor) * erfc(ar) / r );


	}
    }
  }
//final ion, look at the loop indices
t1=NP-1;
  for (int l=0;l<Npol;l++)
    {	
	dpix=pxy[l][0]-xyz[t1][0]; dpix-=dround(dpix/Ls[0])*Ls[0];
	dpiy=pxy[l][1]-xyz[t1][1]; dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));
	//e1(x)=-ei(-x)
	if (rpi2<=rmax2)
   		d3_ci-= Q[t1]*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));
	
	
    }
//cylinder-cylinder as a function of their separation
//should there be a cutoff here too, think about periodicity, will they always be l/2?

if (Npol==2)
	d3_cc=-prefactor*taus[1]*taus[0]*Ls[2]*Exponential_Integral_Ei(-SQR(alpha)*(SQR(pxy[1][0]-pxy[0][0])+SQR(pxy[1][1]-pxy[0][1])));

printf("real d3 %lf d3_ci %lf d3_cc %lf\n", d3,d3_ci,d3_cc);
*e_ii=d3;
*e_pi=d3_ci;
*e_pp=d3_cc;

  return d3+d3_ci+d3_cc;
  //return d3;	
}

/*----------------------------------------------------------------------*/

double Ewald_r_space_line(double** xyz, double *Q, 
		     double *Fx, double *Fy, double *Fz)
{
  int    t1,t2;
  double dx,dy,dz,r,r2,ar,ar2, dpix, dpiy;
  double d1,d2,d3=0.0,d4, d3_ci=0, d3_cc=0;
  double rpi2=0;
  static double wupi = 1.77245385090551602729816748334;

  for (t1 = 0; t1 < NP-1; t1++)   /* Quick and Dirty N^2 loop */
  {
 //cylinder-ion contribution
    for (int l=0;l<Npol;l++)
    {	
	dpix=pxy[l][0]-xyz[t1][0]; dpix-=dround(dpix/Ls[0])*Ls[0];
	dpiy=pxy[l][1]-xyz[t1][1]; dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));
	//e1(x)=-ei(-x)
	if (rpi2<=rmax2)
   		d3_ci-= Q[t1]*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));
	
	
    }
    	

    for (t2 = t1 + 1; t2<NP; t2++) {
      dx = xyz[t1][0] - xyz[t2][0]; dx -= dround(dx/Ls[0])*Ls[0];
      dy = xyz[t1][1] - xyz[t2][1]; dy -= dround(dy/Ls[1])*Ls[1]; 
      dz = xyz[t1][2] - xyz[t2][2]; dz -= dround(dz/Ls[2])*Ls[2];
      r2 = SQR(dx) + SQR(dy) + SQR(dz);
      /*      printf("%d %d %le \n",t1+1,t2+1,sqrt(r2));*/
      /*      printf("%d  %d  %le  %le\n",t1+1,t2+1,r2,rmax2);*/

	if (r2 <= rmax2) {
	  //printf("Forces");
	  ar2 = SQR(ar = alpha * (r = sqrt(r2)));
	  d3 += ( d2 = (d1 = Q[t1] * Q[t2] * prefactor) * erfc(ar) / r );

	  if (Fx||Fy||Fz) {   /* if forces should be computed... */
	    d4  = ( d2 + d1*2.0*alpha*exp(-ar2)/wupi) / r2;
	    
	    Fx[t1] += d4*dx;
	    Fy[t1] += d4*dy;
	    Fz[t1] += d4*dz;
	    Fx[t2] -= d4*dx;
	    Fy[t2] -= d4*dy;
	    Fz[t2] -= d4*dz;
	  }
	}
    }
  }
//final ion, look at the loop indices
t1=NP-1;
  for (int l=0;l<Npol;l++)
    {	
	dpix=pxy[l][0]-xyz[t1][0]; dpix-=dround(dpix/Ls[0])*Ls[0];
	dpiy=pxy[l][1]-xyz[t1][1]; dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));
	//e1(x)=-ei(-x)
	if (rpi2<=rmax2)
   		d3_ci-= Q[t1]*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));
	
	
    }
//cylinder-cylinder as a function of their separation
//should there be a cutoff here too, think about periodicity, will they always be l/2?

if (Npol==2)
	d3_cc=-prefactor*taus[1]*taus[0]*Ls[2]*Exponential_Integral_Ei(-SQR(alpha)*(SQR(pxy[1][0]-pxy[0][0])+SQR(pxy[1][1]-pxy[0][1])));

printf("real line d3 %lf d3_ci %lf d3_cc %lf\n", d3,d3_ci,d3_cc);
  return d3+d3_ci+d3_cc;
  //return d3;	
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Functions for one particle move energy charge by Maria

double Ewald_r_space_part_line(double** xyz, double *Q, 
		     double *MCcoord, double MCQ, int MCion)
{
  int    t2;
  double dx,dy,dz,r,r2,ar, dpix, dpiy, rpi2;
  double d1,d2,d3=0.0;
  double d3_ci=0;

   //cylinder-ion contribution
    for (int l=0;l<Npol;l++)
    {	
	dpix=pxy[l][0]-MCcoord[0]; 
	dpix-=dround(dpix/Ls[0])*Ls[0];
	
	dpiy=pxy[l][1]-MCcoord[1]; 
	dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));

	if (rpi2<=rmax2)
   		d3_ci-= MCQ*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));
	
    }	

    for (t2 = 0; t2<NP; t2++) {
      if ( t2 != MCion )
	{
	  dx = MCcoord[0] - xyz[t2][0]; dx -= dround(dx/Ls[0])*Ls[0];
	  dy = MCcoord[1] - xyz[t2][1]; dy -= dround(dy/Ls[1])*Ls[1]; 
	  dz = MCcoord[2] - xyz[t2][2]; dz -= dround(dz/Ls[2])*Ls[2];
	  r2 = SQR(dx) + SQR(dy) + SQR(dz);

	  
	  if (r2 <= rmax2) {

	    ar = alpha * (r = sqrt(r2));
	    d3 += ( d2 = (d1 = MCQ * Q[t2] * prefactor) * erfc(ar) / r );
	    
	  }
	}
    }

  return d3+d3_ci;
 //  return d3;	
}

/*----------------------------------------------------------------------*/

double Ewald_k_space_part_line(double** xyz, double *Q, double* MCcoord, double MCQ, int MCion)
{
  int    i;
  int    nx,ny,nz;
  double h,c,s,re;
  double d1,d2=0.0,d3=0.0, ci=0.0;
  double hpol[2];	 

  static double PI   = 3.14159265358979323846264338328;
  static double wupi = 1.77245385090551602729816748334;
   
  d1=2*PI/L3;



  sum_q_2 = sum_q2 = 0.0;
  for (i=0; i<NP; i++) {
    sum_q_2 += Q[i];
    sum_q2  += SQR(Q[i]);
  }
  sum_q_2 *= sum_q_2;


	

  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++)
      for (nz=-kmax[2]; nz<=kmax[2]; nz++) {
	//i = MCion;
	e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = 0.0;
  		hpol[0]=2.0 * PI * (nx*pxy[0][0]/Ls[0] + ny*pxy[0][1]/Ls[1]);
 	 //if (Npol==2) removed if for vectorization, pxy sould be initialized
		hpol[1]=2.0 * PI * (nx*pxy[1][0]/Ls[0] + ny*pxy[1][1]/Ls[1]);

	for (i=0; i<NP; i++)
	if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	{
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]) ;
	    e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] += Q[i]*cos(h);
	    e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] -= Q[i]*sin(h);
	
	  }
	if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) {
	  h = 2.0 * PI * (nx*MCcoord[0]/Ls[0] + ny*MCcoord[1]/Ls[1] + nz*MCcoord[2]/Ls[2]) ;
	
	  c = cos(h);
	  s = sin(h);
	  
	  d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
	 // printf("%lf %lf %lf %lf %lf %lf\n", c,s,d2,h,e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]],e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]);	
	  re = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c 
		- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s) * d2;
	
	    if(nz==0)
	    {
		for (int l=0;l<Npol;l++)
			ci=ci+d2*taus[l]*Ls[2]*MCQ*(cos(h)*cos(hpol[l])+sin(h)*sin(hpol[l]));
	    }
         //printf("Q %lf re %lf\n", MCQ, re);
	  d3 += MCQ * re;
	 }
      }

  //printf("%le",d1);
//printf("imaginaariavaruus %lf %lf\n", d3*d1, ci*d1);
  d3=d3+ci;	

	
  d3 *= d1;
  //printf("imaginaariavaruus %lf %lf\n", d3, ci);
  /* self energy and net charge correction: */
  //d3 -= prefactor * ( sum_q2 * alpha / wupi  +  sum_q_2 * PI / (2.0*L3*SQR(alpha)) );
//no net charge correction
  d3 -= prefactor * ( SQR(MCQ) * alpha / wupi );

  /*  printf("\n%le", sum_q2); */
  
  return(2.0*d3);
}
/*-------------------------------------------------------------------------------------*/
//--------------------------------------------------------------------------------------//
//			Optimized ewald functions for two particle move:
//				Quick and ditry implementation
double Ewald_r_space_doublemove(double** xyz, double *Q, 
		     double *MCcoord, double MCQ, int MCion,double *MCcoord2, double MCQ2, int MCion2)
{ 
  int    t2;
  double dx,dy,dz,r,r2,ar, dpix, dpiy, rpi2,rpi2_2, dx2,dy2,dz2,r2_2;
  double d1,d2,d3=0.0;
  double d3_ci=0;

   //cylinder-ion contributions
	
    for (int l=0;l<Npol;l++)
    {	

	dpix=pxy[l][0]-MCcoord[0]; 
	dpix-=dround(dpix/Ls[0])*Ls[0];
	
	dpiy=pxy[l][1]-MCcoord[1]; 
	dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));

	if (rpi2<=rmax2)
   		d3_ci-= MCQ*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2));


	dpix=pxy[l][0]-MCcoord2[0]; 
	dpix-=dround(dpix/Ls[0])*Ls[0];
	
	dpiy=pxy[l][1]-MCcoord2[1]; 
	dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2_2=(SQR(dpix)+SQR(dpiy));

	if (rpi2_2<=rmax2)
   		d3_ci-= MCQ2*prefactor*taus[l]*Exponential_Integral_Ei(-SQR(alpha)*(rpi2_2));
	
    }
    for (t2 = 0; t2<NP; t2++) {
	//mutual contribution of ions dont have to be counted if we only use the difference
      if ( t2 != MCion && t2!=MCion2 )
	{
	  dx = MCcoord[0] - xyz[t2][0]; dx -= dround(dx/Ls[0])*Ls[0];
	  dy = MCcoord[1] - xyz[t2][1]; dy -= dround(dy/Ls[1])*Ls[1]; 
	  dz = MCcoord[2] - xyz[t2][2]; dz -= dround(dz/Ls[2])*Ls[2];
	  r2 = SQR(dx) + SQR(dy) + SQR(dz);

	  dx2 = MCcoord2[0] - xyz[t2][0]; dx2 -= dround(dx2/Ls[0])*Ls[0];
	  dy2 = MCcoord2[1] - xyz[t2][1]; dy2 -= dround(dy2/Ls[1])*Ls[1]; 
	  dz2 = MCcoord2[2] - xyz[t2][2]; dz2 -= dround(dz2/Ls[2])*Ls[2];
	  r2_2 = SQR(dx2) + SQR(dy2) + SQR(dz2);

	  
	  if (r2 <= rmax2) {

	    ar = alpha * (r = sqrt(r2));
	    d3 += ( d2 = (d1 = MCQ * Q[t2] * prefactor) * erfc(ar) / r );
	    
	  }
	  if (r2_2 <= rmax2) {

	    ar = alpha * (r = sqrt(r2_2));
	    d3 += ( d2 = (d1 = MCQ2 * Q[t2] * prefactor) * erfc(ar) / r );
	    
	  }
	}
	
    }

  return d3+d3_ci;
 //  return d3;	
}

/*----------------------------------------------------------------------*/

double Ewald_k_space_doublemove(double** xyz, double *Q, double* MCcoord, double MCQ, double* MCcoord2, double MCQ2)
{
  int    i;
  int    nx,ny,nz;
  double h,c,s,re,h2,c2,s2,re2;
  double d1,d2=0.0,d3=0.0, ci=0.0;
  double hpol[2];	 

  static double PI   = 3.14159265358979323846264338328;
  static double wupi = 1.77245385090551602729816748334;
   
  d1=2*PI/L3;



  sum_q_2 = sum_q2 = 0.0;
  for (i=0; i<NP; i++) {
    sum_q_2 += Q[i];
    sum_q2  += SQR(Q[i]);
  }
  sum_q_2 *= sum_q_2;


	

  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++)
      for (nz=-kmax[2]; nz<=kmax[2]; nz++) {
	//i = MCion;
	e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] = 0.0;
  		hpol[0]=2.0 * PI * (nx*pxy[0][0]/Ls[0] + ny*pxy[0][1]/Ls[1]);
 	 //if (Npol==2) removed if for vectorization, pxy sould be initialized
		hpol[1]=2.0 * PI * (nx*pxy[1][0]/Ls[0] + ny*pxy[1][1]/Ls[1]);


	for (i=0; i<NP; i++)
	if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
	{
	    h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]) ;
	    e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] += Q[i]*cos(h);
	    e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]] -= Q[i]*sin(h);
	
	  }
	if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) {
	  h = 2.0 * PI * (nx*MCcoord[0]/Ls[0] + ny*MCcoord[1]/Ls[1] + nz*MCcoord[2]/Ls[2]) ;
	  h2 = 2.0 * PI * (nx*MCcoord2[0]/Ls[0] + ny*MCcoord2[1]/Ls[1] + nz*MCcoord2[2]/Ls[2]) ;
	  c = cos(h);
	  s = sin(h);
	  c2 = cos(h2);
	  s2 = sin(h2);
	  
	  d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
	 // printf("%lf %lf %lf %lf %lf %lf\n", c,s,d2,h,e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]],e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]);	
	  re = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c 
		- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s) * d2;
	  re2 = (e_re[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*c2 
		- e_im[nx+kmax[0]][ny+kmax[1]][nz+kmax[2]]*s2) * d2;
	
	    if(nz==0)
	    {
		for (int l=0;l<Npol;l++)
		{
			ci=ci+d2*taus[l]*Ls[2]*MCQ*(cos(h)*cos(hpol[l])+sin(h)*sin(hpol[l]));
			ci=ci+d2*taus[l]*Ls[2]*MCQ2*(cos(h2)*cos(hpol[l])+sin(h2)*sin(hpol[l]));
		}	
	    }
         //printf("Q %lf re %lf\n", MCQ, re);
	  d3 += MCQ * re;
	  d3 += MCQ2 * re2;
	 }
      }

  //printf("%le",d1);
  d3=d3+ci;	

	
  d3 *= d1;
  //printf("imaginaariavaruus %lf %lf\n", d3, ci);
  /* self energy and net charge correction: */
  //d3 -= prefactor * ( sum_q2 * alpha / wupi  +  sum_q_2 * PI / (2.0*L3*SQR(alpha)) );
//no net charge correction
  d3 -= prefactor * ( SQR(MCQ) * alpha / wupi );

  /*  printf("\n%le", sum_q2); */
  
  return(2.0*d3);
}

//----------------Force routines for forces on polymers---------------------------------//

void ion_to_cylinder_r(double** xyz,double* Q, int pndx, double* fx, double* fy)
{

double rpi2=0;
double dpix=0;
double dpiy=0;  

for (int t1=0;t1<NP;t1++)
{
	
	dpix=pxy[pndx][0]-xyz[t1][0]; dpix-=dround(dpix/Ls[0])*Ls[0];
	dpiy=pxy[pndx][1]-xyz[t1][1]; dpiy-=dround(dpiy/Ls[1])*Ls[1];
	rpi2=(SQR(dpix)+SQR(dpiy));
	//e1(x)=-ei(-x)
	if (rpi2<=rmax2)
	{
   		*fx+=2*Q[t1]*prefactor*taus[pndx]*exp(-SQR(alpha)*(rpi2))*(dpix/rpi2);
		*fy+=2*Q[t1]*prefactor*taus[pndx]*exp(-SQR(alpha)*(rpi2))*(dpiy/rpi2);
		//*fx+=Q[t1]*prefactor*taus[pndx]*dpix/(rpi2)*(exp(-SQR(alpha)*(rpi2))*2*alpha/sqrt(PI)+erfc(alpha*sqrt(rpi2))/sqrt(rpi2));
		//*fy+=Q[t1]*prefactor*taus[pndx]*dpiy/(rpi2)*(exp(-SQR(alpha)*(rpi2))*2*alpha/sqrt(PI)+erfc(alpha*sqrt(rpi2))/sqrt(rpi2));
	}
	
    

}
}

//----------------------------------------------------------------------------------------
void ion_to_cylinder_k(double** xyz, double* Q,int pndx, double* fx, double* fy)
{

  int    i;
  int    nx,ny,nz;
  double h;
  double d1,d2;
  double hpol;	 

  static double PI   = 3.14159265358979323846264338328;

   
  d1=2*PI/L3;
  nz=0;

  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++) {
	
  		
	for (i=0; i<NP; i++)
	{
		//if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) {

 	  		hpol=2.0 * PI * (nx*pxy[pndx][0]/Ls[0] + ny*pxy[pndx][1]/Ls[1]);
	  		h = 2.0 * PI * (nx*xyz[i][0]/Ls[0] + ny*xyz[i][1]/Ls[1] + nz*xyz[i][2]/Ls[2]);
	  		d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
				
			//*fx+=d1*d2*taus[pndx]*Ls[2]*Q[i]*2*2*PI*nx/Ls[0]*(cos(h)*sin(hpol)-sin(h)*cos(hpol));
			//*fy+=d1*d2*taus[pndx]*Ls[2]*Q[i]*2*2*PI*ny/Ls[1]*(cos(h)*sin(hpol)-sin(h)*cos(hpol));	
			*fx-=d1*d2*taus[pndx]*Ls[2]*Q[i]*2*2*PI*nx/Ls[0]*sin(-hpol+h);
			*fy-=d1*d2*taus[pndx]*Ls[2]*Q[i]*2*2*PI*ny/Ls[1]*sin(-hpol+h);
			
	 	//}
	}
	}

}

//-----------------------------------------------------------------------------------------

void cylinder_to_cylinder_r(double* fx, double* fy, int pndx)
{
//negative direction vector, positive force,. Also positive x-aksis, positive force configurations are generated so that vector from 0 to 1 is negative, ie, pol 0 is at box/2+dist/2 and 1 is at box/2-dist/2. This implementation is not very elegant, but is readable. Could be implemented by using [abs(pndx-1)] 
double rp2=SQR(pxy[0][0]-pxy[1][0])+SQR(pxy[0][1]-pxy[1][1]);
if (rp2<=rmax2)
{
	if (pndx==0)
	{
		*fx=2*prefactor*taus[0]*taus[1]*Ls[2]*exp(-SQR(alpha)*rp2)*((pxy[0][0]-pxy[1][0])/rp2);
		*fy=2*prefactor*taus[0]*taus[1]*Ls[2]*exp(-SQR(alpha)*rp2)*((pxy[0][1]-pxy[1][1])/rp2);	
	}
	if (pndx==1)
	{
		*fx=-2*prefactor*taus[0]*taus[1]*Ls[2]*exp(-SQR(alpha)*rp2)*((pxy[0][0]-pxy[1][0])/rp2);
		*fy=-2*prefactor*taus[0]*taus[1]*Ls[2]*exp(-SQR(alpha)*rp2)*((pxy[0][1]-pxy[1][1])/rp2);		

	}
}

}

//------------------------------------------------------------------------------------------

void cylinder_to_cylinder_k(int pndx, double* fx, double* fy)
{
  
  int    nx,ny,nz;

  double hpol[2];	 
  double d1,d2;
  static double PI   = 3.14159265358979323846264338328;

   
  d1=2*PI/L3;
nz=0;
  for (nx=-kmax[0]; nx<=kmax[0]; nx++)
    for (ny=-kmax[1]; ny<=kmax[1]; ny++) 
	{
		hpol[0]=2.0 * PI * (nx*pxy[0][0]/Ls[0] + ny*pxy[0][1]/Ls[1]);
		hpol[1]=2.0 * PI * (nx*pxy[1][0]/Ls[0] + ny*pxy[1][1]/Ls[1]);	

		 //if (SQR(nx) + SQR(ny) + SQR(nx) < kmax2) 
		//{
	      		d2 = Ghat[abs(nx)][abs(ny)][abs(nz)];
			if (pndx==0)
			{
				*fx+=2*d2*d1*nx*2*PI/Ls[0]*taus[1]*taus[0]*SQR(Ls[2])*sin(hpol[0]-hpol[1]);
				*fy+=2*d2*d1*ny*2*PI/Ls[1]*taus[1]*taus[0]*SQR(Ls[2])*sin(hpol[0]-hpol[1]);
			}
			//vaihda elseen?	
			if (pndx==1)
		 	{
				*fx-=2*d2*d1*nx*2*PI/Ls[0]*taus[1]*taus[0]*SQR(Ls[2])*sin(hpol[0]-hpol[1]);
				*fy-=2*d2*d1*ny*2*PI/Ls[1]*taus[1]*taus[0]*SQR(Ls[2])*sin(hpol[0]-hpol[1]);
			}
		
		//}

	}
}
/*-------------------------------------------------------------*/
//The original ewald function for testing
double Ewald_r_space(double** xyz, double *Q, double** F)
{
  int    t1,t2;
  double dx,dy,dz,r,r2,ar,ar2;
  double d1,d2,d3=0.0,d4;
  double Li=1/Ls[0];
  double Ltmp=Ls[0];  
//printf("L %lf Li %lf pre %lf rmax2 %lf\n", Ltmp,Li, prefactor, rmax2);
  static double wupi = 1.77245385090551602729816748334;

  for (t1 = 0; t1 < NP-1; t1++)   /* Quick and Dirty N^2 loop */
    for (t2 = t1 + 1; t2<NP; t2++) {
      dx = xyz[t1][0] - xyz[t2][0]; dx -= dround(dx*Li)*Ltmp;
      dy = xyz[t1][1] - xyz[t2][1]; dy -= dround(dy*Li)*Ltmp; 
      dz = xyz[t1][2] - xyz[t2][2]; dz -= dround(dz*Li)*Ltmp;
      r2 = SQR(dx) + SQR(dy) + SQR(dz);
      /*      printf("%d %d %le \n",t1+1,t2+1,sqrt(r2));*/
      /*      printf("%d  %d  %le  %le\n",t1+1,t2+1,r2,rmax2);*/

	if (r2 <= rmax2) {
	  //printf("Forces");
	  ar2 = SQR(ar = alpha * (r = sqrt(r2)));
	  d3 += ( d2 = (d1 = Q[t1] * Q[t2] * prefactor) * erfc(ar) / r );
	  printf("alpha %lf t1 %d t2 %d ar2 %lf d3 %lf\n", alpha, t1,t2,ar2,d3);	
	  if (F) {   /* if forces should be computed... */
	    d4  = ( d2 + d1*2.0*alpha*exp(-ar2)/wupi) / r2;
	    
	    F[t1][0] += d4*dx;
	    F[t1][1] += d4*dy;
	    F[t1][2] += d4*dz;
	    F[t2][0] -= d4*dx;
	    F[t2][1] -= d4*dy;
	    F[t2][2] -= d4*dz;
	  }
	}
    }
  return d3;
}
