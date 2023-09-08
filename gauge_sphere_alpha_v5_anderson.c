/*This code is written to study polymer melt confined in a sphere using SCFT. Particular emphasis is on the orientational aspects 
of the chains. Partial saddle-point approximation is used in here. Analytical solution of the divergence free vector order parameter is 
used to construct a vector field in these calculations. 

ONLY  Scalar fields are updated in this code */

/* "Neutral field - N MOD(vector_field)/6" is being updated*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "./util/dnrutil.c"

#include "./util/dnrutil.h"

#include "./util/simpson_int.c"
#include "./util/random_nu.c"
#include "./util/tridag.c"

#define PI 3.1415926535897932384626433832795029
#define tolerance 1e-6

#define NX 32
#define NY 32
#define NZ 32

#define MAXITER 100000
#define NFRAMES 100
#define TOL_MAIN 1.0e-8
#define GAMMA_ADI 0.03          /*Parameter for ADI*/
#define N_ANDER 3

static int RERUN = 0;
static int LIMIT = MAXITER;
static double mult = 5.0;       //Parameter for random number generator used in initial guess
static long iseed = -393;       //seed for random number generator
static int ntotal = NX*NY*NZ;   // Total number of grid points
static int g_total =  NX*NY*NZ; // Total number of guesses
static int n_anchor;            // Total number of guessed points

static double N,chi,RG,Rsphere,para_D,mono_conc,sigma,vol_cell,alpha_para;       /***Parameters for the polymer melt*****/
static double mx_dstep,my_dstep, mz_dstep,ds_step, kl;

static int *indi, NITER;

static double *mx,*my,*mz,***input_adi,***input2_adi,*paramet;
static double *q_s, *qstar_s, *q_s2, *qstar_s2,*logq_qstars,*logq_qstarx, *logq_qstary,*logq_qstarz;
static double *out_f,*x;
static double ***densityP;//,***densityPold,*del_den;
static double ***fieldP,***fPDenom,***poten1,***potential,***logpoten;

/*****New terms ******/
static double ***tx, ***ty, ***tz, ***diva_t, ***curl_tx, ***curl_ty, ***curl_tz;
static double ***fx, ***fy, ***fz;
static double ***psix, ***psiy, ***psiz,***modpsi;
static double Ke1, Ke2, q0, mu;

static double r_x,r_y,r_z;
static double *dA_x,*dA_y,*dA_z,*dB_x,*dB_y,*dB_z,*dC_x,*dC_y,*dC_z;
static double *tempX,*tempY,*tempZ,*dR_x,*dR_y,*dR_z;
static double *d,*dimensions;
static double partition,avgFieldP;
static double **deviations,**fieldP_prev,**Umat, *Vvec;

void initialize() {
  		mx = dvector(1,NX);   /* x */
  		my = dvector(1,NY);   /* y */
  		mz = dvector(1,NZ);   /* z */

  		dimensions = dvector(1,3);  /*Total no of dimensions*/
  		d = dvector(1,2);           /*Dimensions for minimizations*/
  		paramet = dvector(1,2);

  		q_s = dvector(1,ntotal*(NITER+1));
  		qstar_s = dvector(1,ntotal*(NITER+1));
  		q_s2 = dvector(1,ntotal*(NITER+1));
  		qstar_s2 = dvector(1,ntotal*(NITER+1));
  		logq_qstars = dvector(1,ntotal*(NITER+1));
  		logq_qstarx = dvector(1,ntotal*(NITER+1));
  		logq_qstary = dvector(1,ntotal*(NITER+1));
  		logq_qstarz = dvector(1,ntotal*(NITER+1));
    
  		densityP = f3tensor(1,NX,1,NY,1,NZ); 
  		//densityPold = f3tensor(1,NX,1,NY,1,NZ);                
  		//del_den = dvector(1,ntotal);
  		fieldP = f3tensor(1,NX,1,NY,1,NZ);
  		fPDenom = f3tensor(1,NX,1,NY,1,NZ);
  		//fieldS = f3tensor(1,NX,1,NY,1,NZ);          
  		potential = f3tensor(1,NX,1,NY,1,NZ);                  
  		poten1 = f3tensor(1,NX,1,NY,1,NZ);                  
  		logpoten = f3tensor(1,NX,1,NY,1,NZ);                  
  		input_adi = f3tensor(1,NX,1,NY,1,NZ); 
  		input2_adi = f3tensor(1,NX,1,NY,1,NZ); 
  
  		/*****Initialize for all the t's and psis etc *****/
  		tx = f3tensor(1,NX,1,NY,1,NZ);
  		ty = f3tensor(1,NX,1,NY,1,NZ);
  		tz = f3tensor(1,NX,1,NY,1,NZ);
  		fx = f3tensor(1,NX,1,NY,1,NZ);
  		fy = f3tensor(1,NX,1,NY,1,NZ);
  		fz = f3tensor(1,NX,1,NY,1,NZ);
  		diva_t = f3tensor(1,NX,1,NY,1,NZ);
  		curl_tx = f3tensor(1,NX,1,NY,1,NZ);
  		curl_ty = f3tensor(1,NX,1,NY,1,NZ);
  		curl_tz = f3tensor(1,NX,1,NY,1,NZ);
  		psix = f3tensor(1,NX,1,NY,1,NZ);
  		psiy = f3tensor(1,NX,1,NY,1,NZ);
  		psiz = f3tensor(1,NX,1,NY,1,NZ);
  		modpsi = f3tensor(1,NX,1,NY,1,NZ);
  		/***done with all t's and psis etc ******/
              
  		dA_x = dvector(1,NX);
  		dA_y = dvector(1,NY);
  		dA_z = dvector(1,NZ);
              
  		dB_x = dvector(1,NX);
  		dB_y = dvector(1,NY);
  		dB_z = dvector(1,NZ);
              
  		dC_x = dvector(1,NX);
  		dC_y = dvector(1,NY);
  		dC_z = dvector(1,NZ);
              
  		dR_x = dvector(1,NX);
  		tempX = dvector(1,NX);
  		dR_y = dvector(1,NY);
  		tempY = dvector(1,NY);
  		dR_z = dvector(1,NZ);
  		tempZ = dvector(1,NZ);
  
  		indi = ivector(1,1);
  		x = dvector(1,g_total);
  		// x1 = dvector(1,g_total);
  		out_f = dvector(1,g_total);
  		// out_f1 = dvector(1,g_total);
		    
                deviations = dmatrix(1,3,1,ntotal);
		fieldP_prev = dmatrix(1,3,1,ntotal);
		Umat = dmatrix(1,N_ANDER-1,1,N_ANDER-1);
		
		Vvec = dvector(1,N_ANDER-1);
	
        
}

void close_initialize() {
  		free_dvector(mx,1,NX);
  		free_dvector(my,1,NY);
  		free_dvector(mz,1,NZ);
        
  		free_dvector(dimensions,1,3);
  		free_dvector(d,1,2);
 
  		free_dvector(paramet,1,2);
  		free_dvector(q_s,1,ntotal*(NITER+1));
  		free_dvector(qstar_s,1,ntotal*(NITER+1));
  		free_dvector(q_s2,1,ntotal*(NITER+1));
  		free_dvector(qstar_s2,1,ntotal*(NITER+1));
  		free_dvector(logq_qstars,1,ntotal*(NITER+1));
  		free_dvector(logq_qstarx,1,ntotal*(NITER+1));
  		free_dvector(logq_qstary,1,ntotal*(NITER+1));
  		free_dvector(logq_qstarz,1,ntotal*(NITER+1));

  		free_f3tensor(densityP,1,NX,1,NY,1,NZ);
  		// free_f3tensor(densityPold,1,NX,1,NY,1,NZ);
  		// free_dvector(del_den,1,ntotal);  
  		free_f3tensor(fieldP,1,NX,1,NY,1,NZ);
  		free_f3tensor(fPDenom,1,NX,1,NY,1,NZ);
  		//  free_f3tensor(fieldS,1,NX,1,NY,1,NZ);
  		free_f3tensor(potential,1,NX,1,NY,1,NZ);
  		free_f3tensor(poten1,1,NX,1,NY,1,NZ);
  		free_f3tensor(logpoten,1,NX,1,NY,1,NZ);
  		free_f3tensor(input_adi,1,NX,1,NY,1,NZ);
  		free_f3tensor(input2_adi,1,NX,1,NY,1,NZ);

		/******free all t's and psis ********/
		free_f3tensor(tx,1,NX,1,NY,1,NZ);
		free_f3tensor(ty,1,NX,1,NY,1,NZ);
		free_f3tensor(tz,1,NX,1,NY,1,NZ);
  		free_f3tensor(fx,1,NX,1,NY,1,NZ);
		free_f3tensor(fy,1,NX,1,NY,1,NZ);
		free_f3tensor(fz,1,NX,1,NY,1,NZ);
		free_f3tensor(diva_t,1,NX,1,NY,1,NZ);
		free_f3tensor(curl_tx,1,NX,1,NY,1,NZ);	
		free_f3tensor(curl_ty,1,NX,1,NY,1,NZ);
		free_f3tensor(curl_tz,1,NX,1,NY,1,NZ);
		free_f3tensor(psix,1,NX,1,NY,1,NZ);
		free_f3tensor(psiy,1,NX,1,NY,1,NZ);
		free_f3tensor(psiz,1,NX,1,NY,1,NZ);
		free_f3tensor(modpsi,1,NX,1,NY,1,NZ);
		/********freeing done for all ts and psis**************/

  		free_dvector(dA_x,1,NX);
 		free_dvector(dA_y,1,NY);
  		free_dvector(dA_z,1,NZ);
  		free_dvector(dB_x,1,NX);
  		free_dvector(dB_y,1,NY);
  		free_dvector(dB_z,1,NZ);
  		free_dvector(dC_x,1,NX);
  		free_dvector(dC_y,1,NY);
  		free_dvector(dC_z,1,NZ);

  		free_dvector(dR_x,1,NX);
  		free_dvector(tempX,1,NX);
  		free_dvector(dR_y,1,NY);
  		free_dvector(tempY,1,NY);
  		free_dvector(dR_z,1,NZ);
  		free_dvector(tempZ,1,NZ);

  		free_ivector(indi,1,1);

  		free_dvector(x,1,g_total);
		// free_dvector(x1,1,g_total);
  		free_dvector(out_f,1,g_total);

		// free_dvector(out_f1,1,g_total);

		 free_dmatrix(deviations,1,3,1,ntotal);
		 free_dmatrix(fieldP_prev,1,3,1,ntotal);
		 free_dmatrix(Umat,1,N_ANDER-1,1,N_ANDER-1);
		 free_dvector(Vvec,1,N_ANDER-1);

}




void EvolveField(int iTime)
{ int i,j,k,m,t,u,NSTART;
  double dRealT = 0.001;
  double dRealT1 = 0.15;
  double temp,dist;		

                  NSTART = 100;
                 if (iTime <= NSTART)
                   {  
   
                    
                      	for (i=1; i<=NX; i++) {
    			for (j=1; j<=NY; j++) {
      				for (k=1; k<=NZ; k++) {
        					
        						m = k + NZ*((j-1) + NY*(i-1));
							
							out_f[m]  = fieldP[i][j][k] - x[m];
						
                  /*Start saving the values from the last 2 iterations and current one
                       for Anderson mixing*/                   
                 if (iTime == NSTART-2){
                    fieldP_prev[1][m] = fieldP[i][j][k];
                      deviations[1][m] = out_f[m];
                        }
                 else if (iTime == NSTART-1)
                    {
                      fieldP_prev[2][m] = fieldP[i][j][k];
                      deviations[2][m] = out_f[m];
                      }
                    else
                         {
                            fieldP_prev[3][m] = fieldP[i][j][k];
                      deviations[3][m] = out_f[m];

                              }
                  /***********************************/             
                   /*********Simple Mixing**********************/
                   x[m] += dRealT*out_f[m];

		  


                          } /*k*/
			} /*j*/
                         }  /*i*/

                 } /*end of if statement*/
            else
                 {


		    
            /******Computation of V vector and U matrix components************/
                   	/*********U matrix components**************/
                       
			   
		         for(t =1; t<N_ANDER;t++)
                            { 
                           /*******V vectors for Anderson Mixing********/
                                   Vvec[t] = 0.0;
				  
                        for (i=1; i<=NX; i++)
			  {
    			for (j=1; j<=NY; j++)
			  {
      			for (k=1; k<=NZ; k++)
			  {
        					
			  m = k + NZ*((j-1) + NY*(i-1));
							
                        
                                	/*Do the bookkeeping for the next iteration*/
		                 deviations[1][m] = deviations[2][m];
                                 deviations[2][m] = deviations[3][m];
                                 deviations[3][m] = out_f[m];
                                  
                                  fieldP_prev[1][m] = fieldP_prev[2][m];
                                  fieldP_prev[2][m] = fieldP_prev[3][m];
                                  fieldP_prev[3][m] = fieldP[i][j][k];
                                             
                           Vvec[t] +=  (deviations[3][m] - deviations[3-t][m])*deviations[3][m] ;
			   
			    
			   } /*k*/  
                           } /*i*/
                           }  /*j*/
			  	Vvec[t] = Vvec[t]/((double)(ntotal));
				
                          
                     
                
			
			 for(u =1; u<N_ANDER;u++)
                         {
			   Umat[t][u] = 0.0;
                             for (i=1; i<=NX; i++)
			  {
    			for (j=1; j<=NY; j++)
			  {
      			for (k=1; k<=NZ; k++)
			  {
        					
			  m = k + NZ*((j-1) + NY*(i-1));
							
                        
                                	/*Do the bookkeeping for the next iteration*/
		                
			   Umat[t][u] +=  (deviations[3][m] - deviations[3-t][m])*(deviations[3][m] - deviations[3-u][m]); /*For U matrix*/  
			
		
			 } /*k*/  
                           } /*i*/
                           }  /*j*/

			        Umat[t][u] = Umat[t][u]/((double)(ntotal));
                           } /*end of u */
		        
                         } /*end of t */                        

		
			
             /********Calculate Anderson Coefficients********/
                    temp = Umat[1][1]*Umat[2][2]- Umat[2][1]*Umat[1][2];
              Vvec[1] = Vvec[1]*Umat[2][2] - Vvec[2]*Umat[1][2];
              Vvec[1] = Vvec[1]/temp;

             Vvec[2] = Vvec[2]*Umat[1][1] - Vvec[1]*Umat[2][1];
              Vvec[2] = Vvec[2]/temp;

            
               printf("Anderson mixing parameters: %lf\t %lf\n",Vvec[1],Vvec[2]);
          /*dludcmp(UAmat,N_ANDER-1,index,para_lu);
          dlubksb(UAmat,N_ANDER-1,index,VAvec);
            printf("%lf\t %lf\n",VAvec[1],VAvec[2]);
          dludcmp(UBmat,N_ANDER-1,index,para_lu);
          dlubksb(UBmat,N_ANDER-1,index,VBvec);

          dludcmp(USmat,N_ANDER-1,index,para_lu);
          dlubksb(USmat,N_ANDER-1,index,VSvec);*/

                  /*****Guess for next iteration*******/
                   for (i=1; i<=NX; i++)
			  {
    			for (j=1; j<=NY; j++)
			  {
      			for (k=1; k<=NZ; k++)
			  {
        					
			  m = k + NZ*((j-1) + NY*(i-1));
			  
         x[m] = fieldP_prev[3][m] + Vvec[1]*(fieldP_prev[2][m] - fieldP_prev[3][m]); 
        x[m] += Vvec[2]*(fieldP_prev[1][m] - fieldP_prev[3][m]);



 
			  } /*k*/    
 
                            } /*i*/
                           }  /*j*/

                } /*else statement*/

                        }


double single_deri_x(double ***potential, int i, int j, int k) {
		/*Computes one dimensional derivative in x-direction, where value of the quantity (whose derivative is being computed)
		is zero at the boundaries*/
  		double derix;

  		if (i==1) {derix = potential[i+1][j][k];}
  		else if (i==NX) {derix = - potential[i-1][j][k];}
  		else {derix = potential[i+1][j][k] - potential[i-1][j][k];}
  		return derix;
}

double single_deri_y(double ***potential, int i, int j, int k) {
		/*Computes one dimensional derivative in y-direction, where value of the quantity (whose derivative is being computed)
                is zero at the boundaries*/
  		double deriy;

  		if (j==1) {deriy = potential[i][j+1][k];}
  		else if (j==NY) {deriy = - potential[i][j-1][k];}
  		else {deriy = potential[i][j+1][k] - potential[i][j-1][k];}
  		return deriy;
}

double single_deri_z(double ***potential, int i, int j, int k) {
		/*Computes one dimensional derivative in z-direction, where value of the quantity (whose derivative is being computed)
                is zero at the boundaries*/
  		double deriz;
             
  		if (k==1) {deriz = potential[i][j][k+1];}
  		else if (k==NZ) {deriz = - potential[i][j][k-1];}
  		else {deriz = potential[i][j][k+1] - potential[i][j][k-1];}
  		return deriz;
}

double double_deri_x(double ***potential, int i, int j, int k) {
	/*Computes one dimensional Laplacian in x-direction, where value of the quantity (whose derivative is being computed)
                is zero at the boundaries*/
		double derix;
        
  		if (i==1) {derix = potential[i+1][j][k] - 2.0*potential[i][j][k];}
  		else if (i==NX) {derix = -2.0*potential[i][j][k] + potential[i-1][j][k];}
  		else {derix = potential[i+1][j][k] - 2.0*potential[i][j][k] + potential[i-1][j][k];}
  		return derix;
}

double double_deri_y(double ***potential, int i, int j, int k) {
	/*Computes one dimensional Laplacian in y-direction, where value of the quantity (whose derivative is being computed)
                is zero at the boundaries*/
  		double deriy;

  		if (j==1) {deriy = potential[i][j+1][k] - 2.0*potential[i][j][k];}
  		else if (j==NY) {deriy = -2.0*potential[i][j][k] + potential[i][j-1][k];}
  		else {deriy = potential[i][j+1][k] - 2.0*potential[i][j][k] + potential[i][j-1][k];}
  		return deriy;
}

double double_deri_z(double ***potential, int i, int j, int k) {
	/*Computes one dimensional Laplacian in z-direction, where value of the quantity (whose derivative is being computed)
                is zero at the boundaries*/
  		double deriz;
  
  		if (k==1) {deriz = potential[i][j][k+1] - 2.0*potential[i][j][k];}
  		else if (k==NZ) {deriz = -2.0*potential[i][j][k] + potential[i][j][k-1];}
  		else {deriz = potential[i][j][k+1] - 2.0*potential[i][j][k] + potential[i][j][k-1];}
  		return deriz;
}

double triple_deri(double ***potential, int i, int j, int k) {
		/*Computes triple derivatives for the ADI algorithm */
		double deri3;
  		int ap,bp,cp;
  		int an,bn,cn;
  		deri3 = 0.0;

  		/*Take care of dirichlet boundary conditions in all directions*/
  		if (i==NX) {ap = i;} else {ap = i+1;}
  		if (j==NY) {bp = j;} else {bp = j+1;}
  		if (i==1) {an = i;} else {an = i-1;}
 		if (j==1) {bn = j;} else {bn = j-1;}
  		if (k==1) {cn = k;} else {cn = k-1;}
  		if (k==NZ) {cp = k;} else {cp = k+1;}	

  		deri3 += potential[ap][bp][cp] - 2.0*potential[i][bp][cp] + potential[an][bp][cp];
  		deri3 = deri3 - 2.0*(potential[ap][j][cp] - 2.0*potential[i][j][cp] + potential[an][j][cp]);
  		deri3 += potential[ap][bn][cp] - 2.0*potential[i][bn][cp] + potential[an][bn][cp];

  		deri3 = deri3 - 2.0*(potential[ap][bp][k] - 2.0*potential[i][bp][k] + potential[an][bp][k]);
  		deri3 += 4.0*(potential[ap][j][k] - 2.0*potential[i][j][k] + potential[an][j][k]);
  		deri3 = deri3 - 2.0*(potential[ap][bn][k] - 2.0*potential[i][bn][k] + potential[an][bn][k]);

  		deri3 += (potential[ap][bp][cn] - 2.0*potential[i][bp][cn] + potential[an][bp][cn]);
  		deri3 = deri3 - 2.0*(potential[ap][j][cn] - 2.0*potential[i][j][cn] + potential[an][j][cn]);

  		deri3 += (potential[ap][bn][cn] - 2.0*potential[i][bn][cn] + potential[an][bn][cn]);
  
  		return deri3;
}

/********curl and divergence of the vectors represented by ***x ***y ***z ********/
void curlV(double ***ox, double ***oy, double ***oz, double ***ix, double ***iy, double ***iz) {
	int i,j,k;
	for(i=1; i<=NX; i++) {
		for(j=1; j<=NY; j++) {
			for(k=1; k<=NZ; k++) {
			  ox[i][j][k] = single_deri_y(iz,i,j,k)/(2.0*my_dstep) - single_deri_z(iy,i,j,k)/(2.0*mz_dstep);
				oy[i][j][k] = single_deri_z(ix,i,j,k)/(2.0*mz_dstep) - single_deri_x(iz,i,j,k)/(2.0*mx_dstep);
				oz[i][j][k] = single_deri_x(iy,i,j,k)/(2.0*mx_dstep) - single_deri_y(ix,i,j,k)/(2.0*my_dstep);
			}
		}
	}
}

void divV(double ***o, double ***ix, double ***iy, double ***iz) {
	int i,j,k;
	for(i=1; i<=NX; i++) {
		for(j=1; j<=NY; j++) {
			for(k=1; k<=NZ; k++) {
				o[i][j][k] = single_deri_x(ix,i,j,k)/(2.0*mx_dstep) + single_deri_y(iy,i,j,k)/(2.0*my_dstep) + single_deri_z(iz,i,j,k)/(2.0*mz_dstep);
			}
		}
	}	
}

void qlaplacian() {
       /*Computes laplacian of the propagators at each location and chain contour step*/
        int i,j,k,t,m;
	for(t=1; t<=(NITER+1); t++) {
		for(i=1; i<=NX; i++) {
                	for(j=1; j<=NY; j++) {
                        	for(k=1; k<=NZ; k++) { m = k + NZ*((j-1) + NY*(i-1));

							potential[i][j][k] = q_s[(t-1)*ntotal+m];
							poten1[i][j][k] = qstar_s[(NITER-t+1)*ntotal+m];

                        			 if ((potential[i][j][k]==0.0) || (poten1[i][j][k]==0.0)) {logpoten[i][j][k] = 0.0;}
                        			 else {logpoten[i][j][k] = log(fabs(potential[i][j][k])) - log(fabs(poten1[i][j][k]));}
						logq_qstars[(t-1)*ntotal+m] = logpoten[i][j][k] ;
 	
						}
                				}
        				}

                for(i=1; i<=NX; i++) {
                        for(j=1; j<=NY; j++) {
                                for(k=1; k<=NZ; k++) { 
				
					
				//All set for for computing the gradient

				m = k + NZ*((j-1) + NY*(i-1));

                                logq_qstarx[(t-1)*ntotal+m] = single_deri_x(logpoten,i,j,k)/(2.0*mx_dstep);
                                logq_qstary[(t-1)*ntotal+m] = single_deri_y(logpoten,i,j,k)/(2.0*my_dstep);
                                logq_qstarz[(t-1)*ntotal+m] = single_deri_z(logpoten,i,j,k)/(2.0*mz_dstep);

				if((t == 1) || (t == (NITER + 1))) {
          fx[i][j][k] += 0.5*potential[i][j][k]*poten1[i][j][k]*logq_qstarx[(t-1)*ntotal+m]; //single_deri_x(logpoten,i,j,k)/(2.0*mx_dstep); 
          fy[i][j][k] += 0.5*potential[i][j][k]*poten1[i][j][k]*logq_qstary[(t-1)*ntotal+m];//single_deri_y(logpoten,i,j,k)/(2.0*my_dstep); 
          fz[i][j][k] += 0.5*potential[i][j][k]*poten1[i][j][k]*logq_qstarz[(t-1)*ntotal+m]; //single_deri_z(logpoten,i,j,k)/(2.0*mz_dstep); 
          
	   //fx[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_x(poten1,i,j,k)/(2.0*mx_dstep))-poten1[i][j][k]*(single_deri_x(potential,i,j,k)/(2.0*mx_dstep))); 
          //fy[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_y(poten1,i,j,k)/(2.0*my_dstep))-poten1[i][j][k]*(single_deri_y(potential,i,j,k)/(2.0*my_dstep)));
          //fz[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_z(poten1,i,j,k)/(2.0*mz_dstep))-poten1[i][j][k]*(single_deri_z(potential,i,j,k)/(2.0*mz_dstep)));
          							}
	else {

	  fx[i][j][k] += potential[i][j][k]*poten1[i][j][k]*logq_qstarx[(t-1)*ntotal+m]; //single_deri_x(logpoten,i,j,k)/(2.0*mx_dstep);
          fy[i][j][k] += potential[i][j][k]*poten1[i][j][k]*logq_qstary[(t-1)*ntotal+m]; //single_deri_y(logpoten,i,j,k)/(2.0*my_dstep);
          fz[i][j][k] += potential[i][j][k]*poten1[i][j][k]*logq_qstarz[(t-1)*ntotal+m];//single_deri_z(logpoten,i,j,k)/(2.0*mz_dstep);

            //fx[i][j][k] += (potential[i][j][k]*(single_deri_x(poten1,i,j,k)/(2.0*mx_dstep))-poten1[i][j][k]*(single_deri_x(potential,i,j,k)/(2.0*mx_dstep)));
            //fy[i][j][k] += (potential[i][j][k]*(single_deri_y(poten1,i,j,k)/(2.0*my_dstep))-poten1[i][j][k]*(single_deri_y(potential,i,j,k)/(2.0*my_dstep)));
            //fz[i][j][k] += (potential[i][j][k]*(single_deri_z(poten1,i,j,k)/(2.0*mz_dstep))-poten1[i][j][k]*(single_deri_z(potential,i,j,k)/(2.0*mz_dstep)));
          } //else ends 

		
				/*Calculation of laplacians*/
                                
				q_s2[(t-1)*ntotal+m] = double_deri_x(potential,i,j,k)/(mx_dstep*mx_dstep) + double_deri_y(potential,i,j,k)/(my_dstep*my_dstep) + double_deri_z(potential,i,j,k)/(mz_dstep*mz_dstep);
                                qstar_s2[(t-1)*ntotal+m] = double_deri_x(poten1,i,j,k)/(mx_dstep*mx_dstep) + double_deri_y(poten1,i,j,k)/(my_dstep*my_dstep) + double_deri_z(poten1,i,j,k)/(mz_dstep*mz_dstep);
                                                }
                                                }
                                        }


			}/*t*/
}

/**********curl and divergence of the vectors end here **************/

void StartArrage_x() {
  		int i;
  		for(i=1; i<=NX; i++) {
   					dA_x[i] = -0.5*(1.0 + GAMMA_ADI)*r_x;
					dB_x[i] = 1.0 + (1.0 + GAMMA_ADI)*r_x;
    					dC_x[i] = dA_x[i];
    					// printf("x %d\t %e\t %e\t %e\n",i,dA_x[i],dB_x[i],dC_x[i]);
					}
}

void StartArrage_y() { 
  		int i;

		for(i=1; i<=NY; i++) {
   					dA_y[i] = -0.5*(1.0 + GAMMA_ADI)*r_y;
					dB_y[i] = 1.0 + (1.0 + GAMMA_ADI)*r_y;
    					dC_y[i] = dA_y[i];
    // printf("y %d\t %e\t %e\t %e\n",i,dA_y[i],dB_y[i],dC_y[i]);
 					}
}

void StartArrage_z() { 
  		int i;

		for(i=1; i<=NZ; i++) {
    					dA_z[i] = -0.5*(1.0 + GAMMA_ADI)*r_z;
					dB_z[i] = 1.0 + (1.0 + GAMMA_ADI)*r_z;
    					dC_z[i] = dA_z[i];
    //printf("z %d\t %e\t %e\t %e\n",i,dA_z[i],dB_z[i],dC_z[i]);
					}
}

void UpArrage_x(double ***potential, double ***input_adi,int j, int k) { 
		int i,m;
  		double trip_d;
  
  		for(i=1; i<=NX; i++) {
   					dR_x[i] = potential[i][j][k];

 					dR_x[i] += (1.0 - 0.5*(1.0 + GAMMA_ADI))*r_x*double_deri_x(potential,i,j,k);
 					dR_x[i] += r_y*double_deri_y(potential,i,j,k);
 					dR_x[i] += r_z*double_deri_z(potential,i,j,k);
 
					trip_d = triple_deri(potential,i,j,k);
					dR_x[i] += (6.0*GAMMA_ADI + 2.0*pow(GAMMA_ADI,3.0))*r_x*r_y*r_z*trip_d/8.0;
    					// m = k + NZ*((j-1) + NY*(i-1));

					dR_x[i] += input_adi[i][j][k];
 					}
}

void UpArrage_y(double ***potential, double ***poten1,int i, int k) { 
		int j;
  		for(j=1; j<=NY; j++) {
   					dR_y[j] = poten1[i][j][k]; /*field  at n + 1/3*/
    					dR_y[j] = dR_y[j] - 0.5*(1.0 + GAMMA_ADI)*r_y*double_deri_y(potential,i,j,k);
  					}
}

void UpArrage_z(double ***potential, double ***poten1,int i, int j) { 
  		int k;
  		for(k=1; k<=NZ; k++) {
   					dR_z[k] = poten1[i][j][k]; /*potential at n + 2/3*/
					dR_z[k] = dR_z[k] - 0.5*(1.0 + GAMMA_ADI)*r_z*double_deri_z(potential,i,j,k);
					}
}


void solver_adi( double ***potential, double ***input_adi) {
		/*Solves for propagator at the next time step starting from charge density like term given by input_adi
 		The routine uses 3D Douglas ADi with gamma parameter.*/
  		int i,j,k,m;

  		/*ADI - First step along X*/ 
  		for (j=1; j<=NY; j++) {
    			for (k=1; k<=NZ; k++) {
      						UpArrage_x(potential,input_adi,j,k); 
      		/*Solve the tridiagonal set of equations arising because of density 
		  being zero at the boundaries*/ 
      						tridag(dA_x,dB_x,dC_x,dR_x,tempX,NX);      
      		/*Save the results at 1/3 of the step in another variable*/  
      		for (i=1; i<=NX; i++) {
        				poten1[i][j][k] = tempX[i]; 
        	//printf("po1 %d\t %d\t %d\t %e\t %e\t %e\n",i,j,k,poten1[i][j][k],potential[i][j][k],input_adi[i][j][k]) ;
      					}
    						}
  					}                   

  		/*ADI - Second step along Y*/
  		for (i=1; i<=NX; i++) {
    			for (k=1; k<=NZ; k++) {
      						UpArrage_y(potential,poten1,i,k);
      		/*Solve the tridiagonal set of equations arising because of density 
		  being zero at the boundaries*/ 
      						tridag(dA_y,dB_y,dC_y,dR_y,tempY,NY);                           
      		/*Save the results at 2/3 of the step in the variable used to store 
		  results at 1/3 step - this is possible because derivative of 1/3 variables are 
		  not involved in the calculation*/  
      		for (j=1; j<=NY; j++) {
        				poten1[i][j][k] = tempY[j]; /*n+2/3*/
      					}
    						}
  					}

  		/*ADI - Third step along Z*/
  		for (i=1; i<=NX; i++) {
   	 		for (j=1; j<=NY; j++) {
      						UpArrage_z(potential,poten1,i,j);
      		/*Solve the tridiagonal set of equations arising because of density 
		  being zero at the boundaries*/ 
      						tridag(dA_z,dB_z,dC_z,dR_z,tempZ,NZ);     
      			for (k=1; k<=NZ; k++) {
        					potential[i][j][k] = tempZ[k]; /*n+1 i.e., potential[i][j][k] = potentialnew[i][j][k]; */
        					//printf("%d\t %d\t %d\t %lf\n",i,j,k,potential[i][j][k]);
      						}
    						}
  					}         
  		return;
}

double single_deri_prop(double *prop, int t, int m) {
                /*First derivative of propagator in time. Not periodic w.r.t. chain contour*/
                double deri_prop;
                if (t==0) {deri_prop = 0.0;} /*Tension is taken to be zero at the ends*/
                else if (t==NITER) {deri_prop = 0.0;} /*Tension is taken to be zero at the ends*/
                else {deri_prop = prop[(t+1)*ntotal+m]-prop[(t-1)*ntotal+m];} /*Central difference in time*/
                return deri_prop;
}


int CalDensity_FieldP(double ***densityP,double ***fieldP,double *out_f,double *x) { 
      /*Calculates the monomer density(divided by the density if the monomers are distributed uniformly), new field and difference between old and new fields */
		int i,j,k,t,m;
  		int iStable = 1;
		double tmp1;
		double dRealT = 0.001;
                double dRealT1 = 0.0;
                double dDelta1 = 0.0;
                double dDelta2 = 0.0;
                double dDmax = 0.0;
		/*Compute laplacian of the propagators at each time step and grid point*/
                qlaplacian();

		partition = 0.0;
	//	avgFieldP = 0.0;
  		for (i=1; i<=NX; i++) {
    			for (j=1; j<=NY; j++) {
      				for (k=1; k<=NZ; k++) {
        						/*Calculation of monomer density*/
        						m = k + NZ*((j-1) + NY*(i-1));
        						densityP[i][j][k] = 0.0;
		    					partition += mx_dstep*my_dstep*mz_dstep*q_s[(NITER)*ntotal+m]; 
							//fieldP[i][j][k] = qstar_s[(NITER)*ntotal+m] - q_s[(NITER)*ntotal+m];
							fieldP[i][j][k] = 0.0;
                                                        fPDenom[i][j][k] = 0.0;
                                                        densityP[i][j][k] = 0.0;

        						/*Integration over chain contour variable - t*/
        						for(t=1; t<=(NITER+1); t++) { 
									tmp1 = q_s[(t-1)*ntotal+m]*qstar_s[(NITER+1-t)*ntotal+m];
          							if((t == 1) || (t == (NITER + 1))) 
									{
									densityP[i][j][k] += 0.5*tmp1;//q_s[(t-1)*ntotal+m]*qstar_s[(NITER+1-t)*ntotal+m];

 				/* This expression of the field results from the partial saddle-point approximation */

	//fieldP[i][j][k] -= 0.5*(2.0*q_s[(t-1)*ntotal+m]*single_deri_prop(qstar_s,NITER+1-t,m))/(2.0*ds_step);
	fieldP[i][j][k] -= 0.5*tmp1*single_deri_prop(logq_qstars,t-1,m);
	//fieldP[i][j][k] -= 0.5*(qstar_s[(NITER+1-t)*ntotal+m]*single_deri_prop(q_s,t-1,m)+q_s[(t-1)*ntotal+m]*single_deri_prop(qstar_s,NITER+1-t,m))/(2.0*ds_step);
                 fieldP[i][j][k] += 0.5*(qstar_s[(NITER+1-t)*ntotal+m]*q_s2[(t-1)*ntotal+m] + q_s[(t-1)*ntotal+m]*qstar_s2[(NITER+1-t)*ntotal+m]);

                 //fPDenom[i][j][k] += 0.5*q_s[(t-1)*ntotal+m]*qstar_s[(NITER+1-t)*ntotal+m];

									}
          							else {
									densityP[i][j][k] += tmp1;//q_s[(t-1)*ntotal+m]*qstar_s[(NITER+1-t)*ntotal+m];
	//fieldP[i][j][k] -= (2.0*q_s[(t-1)*ntotal+m]*single_deri_prop(qstar_s,NITER+1-t,m))/(2.0*ds_step);
	fieldP[i][j][k] -= tmp1*single_deri_prop(logq_qstars,t-1,m); 
	//fieldP[i][j][k] -= (qstar_s[(NITER+1-t)*ntotal+m]*single_deri_prop(q_s,t-1,m)+q_s[(t-1)*ntotal+m]*single_deri_prop(qstar_s,NITER+1-t,m))/(2.0*ds_step);
                fieldP[i][j][k] += (qstar_s[(NITER+1-t)*ntotal+m]*q_s2[(t-1)*ntotal+m] + q_s[(t-1)*ntotal+m]*qstar_s2[(NITER+1-t)*ntotal+m]);
                 //fPDenom[i][j][k] += q_s[(t-1)*ntotal+m]*qstar_s[(NITER+1-t)*ntotal+m];

									}
        										} /*t*/    
      						
						if (densityP[i][j][k]<=1E-6) {fieldP[i][j][k] = 0.0;}
            					else {fieldP[i][j][k] = fieldP[i][j][k]/(2.0*densityP[i][j][k]);
						
			fieldP[i][j][k] -= (psix[i][j][k]*fx[i][j][k]+psiy[i][j][k]*fy[i][j][k]+psiz[i][j][k]*fz[i][j][k])/densityP[i][j][k];	
							}
			 /*Compute differences between guesses and new quantities*/
                                //m = k + NZ*((j-1) + NY*(i-1));
                                out_f[m] = fieldP[i][j][k] - x[m];
				// x[m] = x[m] + dRealT*out_f[m];

                                dDelta1 = fabs(out_f[m]);
                                if(dDelta1>dDmax) {dDmax = dDelta1;}


							}/*k*/
    						} /*j*/
  					} /*i*/

	    	
		printf("%lf\t",dDmax);
                if (dDmax>tolerance) {iStable = 0;}
                else {iStable=1;}

                printf("%d\n",iStable);
                return iStable;

} 

double ubar(double r, double lambda)
{       /*u(r,lambda,theta) = ubar(r,lambda)*cos(theta) and this routine evaluates ubar(r,lambda)*/

		double ubarfunc;
        	ubarfunc = sin(lambda*r)/(lambda*r);
        	ubarfunc = ubarfunc -cos(lambda*r);
        	ubarfunc = 2.0*ubarfunc/pow(lambda*r,2.0);
		return ubarfunc; 
}

double vbar(double r, double lambda)
{       /*v(r,lambda,theta) = vbar(r,lambda)*sin(theta) and this routine evaluates vbar(r,lambda)*/

        	double vbarfunc;
        	vbarfunc = sin(lambda*r)/(lambda*r);
        	vbarfunc = cos(lambda*r) - vbarfunc;
        	vbarfunc += lambda*r*sin(lambda*r);
        	vbarfunc =  -1.0*vbarfunc/(lambda*r);
        	return vbarfunc;
}

double wbar(double r, double lambda)
{       /*w(r,lambda,theta) = wbar(r,lambda)*sin(theta) and this routine evaluates wbar(r,lambda)*/

        	double wbarfunc;
        	wbarfunc = sin(lambda*r)/(lambda*r);
        	wbarfunc = wbarfunc - cos(lambda*r);
        	wbarfunc =  wbarfunc/(lambda*r);
        	return wbarfunc;
}

double talong_x(double x, double y, double z, double lambda) {
  	/*x-component of a vector transformed from spherical polar to cartesian. The vector in spherical polar co-ordinates is of the form: 

	  t(r,theta,phi,lambda) = ubar(r,lambda)*cos(theta) rhat + vbar(r,lambda)*sin(theta) thetahat + wbar(r,lambda)*sin(theta) phihat
	  rhat, thetahat and phihat are unit vector along r,theta and phi, respectively. 	
		*/
	
		double txcomp,numer,radial;
		radial = sqrt(x*x + y*y + z*z);

		if (radial==0.0) {return 0.0;}
		else 	
		{
    		txcomp = ubar(radial,lambda) + vbar(radial,lambda);
		txcomp = txcomp*x*z/(radial*radial);
                txcomp =  txcomp - wbar(radial,lambda)*(y/radial);
	  	return txcomp;
		}
}

double talong_y(double x, double y, double z, double lambda) {
	/*y-component of a vector transformed from spherical polar to cartesian. The vector in spherical polar co-ordinates is of the form: 

          t(r,theta,phi,lambda) = ubar(r,lambda)*cos(theta) rhat + vbar(r,lambda)*sin(theta) thetahat + wbar(r,lambda)*sin(theta) phihat
          rhat, thetahat and phihat are unit vector along r,theta and phi, respectively.        
                */
		
		double tycomp,numer,radial;
		radial = sqrt(x*x + y*y + z*z);
	
		if (radial==0.0) {return 0.0;}
		else 
		{
	  	tycomp = ubar(radial,lambda) + vbar(radial,lambda);
		tycomp = tycomp*y*z/(radial*radial);
		tycomp += wbar(radial,lambda)*x/radial;
    		return tycomp;

		}
}

double talong_z(double x, double y, double z, double lambda) {
 	/*z-component of a vector transformed from spherical polar to cartesian. The vector in spherical polar co-ordinates is of the form: 

          t(r,theta,phi,lambda) = ubar(r,lambda)*cos(theta) rhat + vbar(r,lambda)*sin(theta) thetahat + wbar(r,lambda)*sin(theta) phihat
          rhat, thetahat and phihat are unit vector along r,theta and phi, respectively.        
                */
	
		double tzcomp,numer,radial;
		radial = sqrt(x*x + y*y + z*z);	
		if (radial==0.0) {return 0.0;}//(2.0/3.0);}
		else 
		{
    		tzcomp = ubar(radial,lambda)*z*z/(radial*radial);
		tzcomp -= vbar(radial,lambda)*x*x/(radial*radial);
		tzcomp -= vbar(radial,lambda)*y*y/(radial*radial);
    		return tzcomp;
			
		}
}

void calc_t_geometry() {

		int i,j,k,t,m,t1;
		for (i=1; i<=NX; i++) {
    			for (j=1; j<=NY; j++) {
      				for (k=1; k<=NZ; k++) {
        						tx[i][j][k] = talong_x(mx[i],my[j],mz[k],alpha_para/Rsphere);
        						ty[i][j][k] = talong_y(mx[i],my[j],mz[k],alpha_para/Rsphere);
        						tz[i][j][k] = talong_z(mx[i],my[j],mz[k],alpha_para/Rsphere);
		  				//printf("%lf\t%lf\t%lf\t%lf\t %lf\t %lf\n",mx[i],my[j],mz[k],tx[i][j][k],ty[i][j][k],tz[i][j][k]);
							}
						}
  					}
}

void calc_f_scft() {
  int i,j,k,t,m,t1;
  double mag=0;

  //gradient of q_s is found using potential and gradient of qstar_s is found using poten1
  for (t=1; t<=(NITER+1); t++) {
    for (i=1; i<=NX; i++) {
      for (j=1; j<=NY; j++) {
        for (k=1; k<=NZ; k++) {
          m = k + NZ*((j-1) + NY*(i-1));
          potential[i][j][k] = q_s[(t-1)*ntotal+m];
          poten1[i][j][k] = qstar_s[(t-1)*ntotal+m];
        }
      }
    }
    //All set for for computing the gradient
        
    for (i=1; i<=NX; i++) {
      for (j=1; j<=NY; j++) {
        for (k=1; k<=NZ; k++) {
          if((t == 1) || (t == (NITER + 1))) {
            fx[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_x(poten1,i,j,k)/(2.0*mx_dstep))-poten1[i][j][k]*(single_deri_x(potential,i,j,k)/(2.0*mx_dstep))); 
            fy[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_y(poten1,i,j,k)/(2.0*my_dstep))-poten1[i][j][k]*(single_deri_y(potential,i,j,k)/(2.0*my_dstep)));
            fz[i][j][k] += 0.5*(potential[i][j][k]*(single_deri_z(poten1,i,j,k)/(2.0*mz_dstep))-poten1[i][j][k]*(single_deri_z(potential,i,j,k)/(2.0*mz_dstep)));
          }
          else {
            fx[i][j][k] += (potential[i][j][k]*(single_deri_x(poten1,i,j,k)/(2.0*mx_dstep))-poten1[i][j][k]*(single_deri_x(potential,i,j,k)/(2.0*mx_dstep)));
            fy[i][j][k] += (potential[i][j][k]*(single_deri_y(poten1,i,j,k)/(2.0*my_dstep))-poten1[i][j][k]*(single_deri_y(potential,i,j,k)/(2.0*my_dstep)));
            fz[i][j][k] += (potential[i][j][k]*(single_deri_z(poten1,i,j,k)/(2.0*mz_dstep))-poten1[i][j][k]*(single_deri_z(potential,i,j,k)/(2.0*mz_dstep))); 
          }
          //else ends 
        }
      }
    }
  }
  
  //f times R_g^3 divided by rho_h times R_g^3, where \rho_h = n N/V 
  for (i=1; i<=NX; i++) {
    for (j=1; j<=NY; j++) {
      for (k=1+n_anchor; k<=NZ; k++) {
        fx[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(6.0*partition);
        fy[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(6.0*partition);
        fz[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(6.0*partition);
		  }
		}
  }	
}

void calc_psi() {
    	//Calculate the vector potential (= -3*t_vec*sqrt(N) ) and the multiplicative factor of sqrt(N) arises in non-dimensional form 
	//of the equations for the propagators
	int k,i,j;
  	double prefac;
	prefac =  Ke2*(alpha_para/Rsphere)*(alpha_para/Rsphere) + 2.0*Ke2*q0*(alpha_para/Rsphere) - mu;
	
	FILE *chiral_field_final1;
        char fname_chiral_field_final1[40];
        sprintf(fname_chiral_field_final1,"chiral_field_final.dat");
        chiral_field_final1= fopen(fname_chiral_field_final1,"w");

        FILE *chiral_vector_final1;
        char fname_chiral_vector_final1[40];
        sprintf(fname_chiral_vector_final1,"chiral_vector_final.dat");
        chiral_vector_final1=fopen(fname_chiral_vector_final1,"w");

	
	for (i=1; i<=NX; i++) {
        	for (j=1; j<=NY; j++) {
      			for (k=1; k<=NZ; k++) {
        					psix[i][j][k] = sqrt(N)*prefac*tx[i][j][k];
        					psiy[i][j][k] = sqrt(N)*prefac*ty[i][j][k];
        					psiz[i][j][k] = sqrt(N)*prefac*tz[i][j][k];
		  				
        					//psix[i][j][k] = -3.0*sqrt(N)*prefac*tx[i][j][k];
        					//psiy[i][j][k] = -3.0*sqrt(N)*prefac*ty[i][j][k];
        					//psiz[i][j][k] = -3.0*sqrt(N)*prefac*tz[i][j][k];
						modpsi[i][j][k] = (psix[i][j][k]*psix[i][j][k]);
						modpsi[i][j][k]  += (psiy[i][j][k]*psiy[i][j][k]);
						modpsi[i][j][k]  += (psiz[i][j][k]*psiz[i][j][k]); 
		  			//	printf("%lf\t%lf\t%lf\t%lf\n",mx[i],my[j],mz[k],modpsi[i][j][k]);
						
	fprintf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], psix[i][j][k],psiy[i][j][k],psiz[i][j][k]);
	//fprintf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], psix[i][j][k],psiy[i][j][k],psiz[i][j][k],modpsi[i][j][k]);
        fprintf(chiral_vector_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], tx[i][j][k],ty[i][j][k],tz[i][j][k]);

						}
					}
  				}	                      

                fclose(chiral_field_final1);
                fclose(chiral_vector_final1);	

}

void calc_psi_scft() {
  int k,i,j;

  divV(diva_t,tx,ty,tz);
  curlV(curl_tx, curl_ty, curl_tz,tx,ty,tz);

  for (i=1; i<=NX; i++) {
    for (j=1; j<=NY; j++) {
      for (k=1; k<=NZ; k++) {
        //grad term
        psix[i][j][k] =- Ke1*single_deri_x(diva_t,i,j,k)/(2.0*mx_dstep); //Ke=K1/RG^2
        psiy[i][j][k] =- Ke1*single_deri_y(diva_t,i,j,k)/(2.0*my_dstep);
        psiz[i][j][k] =- Ke1*single_deri_z(diva_t,i,j,k)/(2.0*mz_dstep);

        //Curl term
        psix[i][j][k] += -Ke2*single_deri_y(curl_tz,i,j,k)/(2.0*my_dstep)+Ke2*single_deri_z(curl_ty,i,j,k)/(2.0*mz_dstep);
        psiy[i][j][k] += -Ke2*single_deri_z(curl_tx,i,j,k)/(2.0*mz_dstep)+Ke2*single_deri_x(curl_tz,i,j,k)/(2.0*mx_dstep);
        psiz[i][j][k] += -Ke2*single_deri_x(curl_ty,i,j,k)/(2.0*mx_dstep)+Ke2*single_deri_y(curl_tx,i,j,k)/(2.0*my_dstep);
  
		    //helix term
		    psix[i][j][k] += 2.0*Ke2*q0*curl_tx[i][j][k];
        psiy[i][j][k] += 2.0*Ke2*q0*curl_ty[i][j][k];
        psiz[i][j][k] += 2.0*Ke2*q0*curl_tz[i][j][k];
      }
    }
  }
}

double output_free() {
  		/*Free energy*/
		int i,j,k;
  		double heint1=0.0;
  		double heint2=0.0;
  		double feint1=0.0;
  		double feint2=0.0;
  		double feterm=0.0;
  		double prefac;
        	prefac =  Ke2*(alpha_para/Rsphere)*(alpha_para/Rsphere) + 2.0*Ke2*q0*(alpha_para/Rsphere) - mu;

		for (i=1; i<=NX; i++) {
    			for (j=1; j<=NY; j++) {
      				for (k=1; k<=NZ; k++) {
				 		/*Here, we are computing rho_p(r)/rho_p0, and rho_p0 = nN/V*/
                                                densityP[i][j][k] *= ds_step*4.0*M_PI*dimensions[1]*dimensions[2]*dimensions[3]/(3.0*partition);
                                                        //printf("%lf\t%lf\t%lf\t%lf\n",mx[i],my[j],mz[k],densityP[i][j][k]);
                                                //fx[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(sqrt(6.0)*partition);
                                                //fy[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(sqrt(6.0)*partition);
                                                //fz[i][j][k]*=ds_step*dimensions[1]*dimensions[2]*dimensions[3]/(sqrt(6.0)*partition);

        				heint1 += mx_dstep*my_dstep*mz_dstep*densityP[i][j][k]*densityP[i][j][k];
        		heint2 += mx_dstep*my_dstep*mz_dstep*N*(tx[i][j][k]*tx[i][j][k] + ty[i][j][k]*ty[i][j][k] + tz[i][j][k]*tz[i][j][k]);
		    	feint1 += mx_dstep*my_dstep*mz_dstep*(fieldP[i][j][k] + (modpsi[i][j][k]/6.0))*densityP[i][j][k]; 
		    	feint2 += mx_dstep*my_dstep*mz_dstep*prefac*sqrt(N)*(tx[i][j][k]*psix[i][j][k] + ty[i][j][k]*psiy[i][j][k] + tz[i][j][k]*psiz[i][j][k]);

      }
    }
  }
	heint1 = 3.0*heint1/(4.0*M_PI*dimensions[1]*dimensions[2]*dimensions[3]);
	heint2 = 3.0*heint2/(4.0*M_PI*dimensions[1]*dimensions[2]*dimensions[3]);
	feint1 = 3.0*feint1/(4.0*M_PI*dimensions[1]*dimensions[2]*dimensions[3]);
	feint2 = 0.5*3.0*feint2/(4.0*M_PI*dimensions[1]*dimensions[2]*dimensions[3]);
     
  	if (partition<=0.0){feterm = 1.0e+7;}
	else{feterm = log(3.0*partition/(4.0*M_PI*dimensions[2]*dimensions[2]*dimensions[3]));}

 
  printf("\n int (rho_p^2) = %lf\n",heint1);
  printf("(int (t^2) = %lf\n",heint2);
  printf("-int(w_p*rho_p) = %lf\n",-feint1);
  printf("-0.5*int(W_p*t_p) = %lf\n",-feint2);
  printf("-log(Q_p) = %lf\n",-feterm);
  printf("Q_p = %lf\n",partition);
  printf("TOTAL (H-H0)/nkT = %lf\n",-feint1 - feint2 - feterm);
  return (- feint1 - feint2 - feterm);

}

int vecfun(int n, double *d, double *x, double *out_f) {
  		
		int k,i,j,m,time,iDone;
  		double totalS,totalC;
  		int ap,an,bp,bn,cp,cn;
  		double prefactor,local_elec;

  		/*******Calculation of q and qstar at all time steps using 
 		 initial conditions and 3D Douglas ADI**********/
            
  		for (time=1; time<=NITER; time++) {
			for(i=1; i<=NX; i++) {
		  		for(j=1; j<=NY; j++) {
		    			for(k=1; k<=NZ; k++) {
		      						m = k + NZ*((j-1) + NY*(i-1) );
          							/*Initialization of q_s to be solved by ADI*/
          							potential[i][j][k] = q_s[(time-1)*ntotal + m];
		    						 /*Initialization of qstar_s to be solved by ADI*/
                        					poten1[i][j][k] = qstar_s[(time-1)*ntotal + m];

								}
      							}
    						}
		/*Prepare to initialize the variables for ADI solver*/
    		for (i=1; i<=NX; i++) {
      			for (j=1; j<=NY; j++) {
        			for (k=1; k<=NZ; k++) {
          						m = k + NZ*((j-1) + NY*(i-1) );
          						/*Charge density like term in ADI solver for q_s*/	
		      					input_adi[i][j][k] = -ds_step*x[m]*q_s[(time-1)*ntotal + m];
          					input_adi[i][j][k] += ds_step*sqrt(2.0/3.0)*(psix[i][j][k]*single_deri_x(potential,i,j,k)/(2.0*mx_dstep)+psiy[i][j][k]*single_deri_y(potential,i,j,k)/(2.0*my_dstep)+psiz[i][j][k]*single_deri_z(potential,i,j,k)/(2.0*mz_dstep)); 
        						

							/*Charge density like term in ADI solver for qstar_s*/
						 /*qstar_s has a negative sign for the vector field*/

          					 input2_adi[i][j][k] = -ds_step*x[m]*qstar_s[(time-1)*ntotal + m];
          					input2_adi[i][j][k] -= ds_step*sqrt(2.0/3.0)*(psix[i][j][k]*single_deri_x(poten1,i,j,k)/(2.0*mx_dstep)+psiy[i][j][k]*single_deri_y(poten1,i,j,k)/(2.0*my_dstep)+psiz[i][j][k]*single_deri_z(poten1,i,j,k)/(2.0*mz_dstep));
						
							}
      						}
    					}

    		/*Call the ADI solver to compute q_s at the next time step*/
    		solver_adi(potential,input_adi);
    		/*Call the ADI solver to compute qstar_s at the next time step*/
    		solver_adi(poten1,input2_adi);
    
    		/*Save the computed quantities so that they can be used for the next time step*/
    		for (i=1; i<=NX; i++) {
      			for (j=1; j<=NY; j++) {
        			for (k=1; k<=NZ; k++) {
          						m = k + NZ*((j-1) + NY*(i-1) );
          						q_s[time*ntotal + m] = potential[i][j][k];
          						qstar_s[time*ntotal + m] = poten1[i][j][k];
         						 /* Spherical boundary condition */
          	if (mx[i]*mx[i]+my[j]*my[j]+mz[k]*mz[k] >= dimensions[3]*dimensions[3]) {q_s[time*ntotal+m] = 0.0;qstar_s[time*ntotal + m] = 0.0;}
        						}
      						}
    					}
  	} /*End of solution for q_s and qstar_s for all times*/

  /*Compute new densities and scalar fields*/
  iDone = CalDensity_FieldP(densityP,fieldP,out_f,x);
  return iDone;
}

double min_func_char(double *d, double *paramet) {
		/*This routine returns the free energy of the system*/ 
		int i,j,k,m,p;
		int ICOUNT;
		int iDONE;
		double free_energy, temp, mag, bw, hstar;

    			para_D = paramet[1];
    			sigma = paramet[2];
    
    			dimensions[1] = d[1];       
    			dimensions[2] = d[2];       
    			vol_cell = dimensions[1]*dimensions[2]*dimensions[3];     

    			/* Set up the grid and grid parameters */
    			mx[1] = -dimensions[1];
    			my[1] = -dimensions[2];
    			mz[1] = -dimensions[3];
    			mx_dstep = 2.0*dimensions[1]/((double)(NX-1));
    			my_dstep = 2.0*dimensions[2]/((double)(NY-1));
    			mz_dstep = 2.0*dimensions[3]/((double)(NZ-1));
       
    		for (p=2; p<=NX; p++) {
      					mx[p] = mx[p-1] + mx_dstep;
    					}
    		for (p=2; p<=NY; p++) {
      					my[p] = my[p-1] + my_dstep;
    					}
    		for (p=2; p<=NZ; p++) {
      					mz[p] = mz[p-1] + mz_dstep;
    					}

    		/****Set up parameters for 3D Douglas ADI*******/
    		r_x = ds_step/pow(mx_dstep,2.0); 
    		r_y = ds_step/pow(my_dstep,2.0); 
    		r_z = ds_step/pow(mz_dstep,2.0); 
    		//printf("%lf\n %lf\n %lf\n",r_x,r_y,r_z);
    
    		StartArrage_x();
    		StartArrage_y();
    		StartArrage_z();

    		/**********************************/
    		/* Input for Non Linear Solver */
	  	k = 64;
	  	printf("limit height=%lf\n",(1.0*k*dimensions[3]/(1.0*NZ)));

	  	/******* Vector field for spherical geometry: psi *******/
    		calc_t_geometry();
    		calc_psi();
      
		FILE *chiral_field_final1;
        	char fname_chiral_field_final1[1000];
        	
		if (RERUN == 0){sprintf(fname_chiral_field_final1,"chiral_field_final.dat");}
		else{sprintf(fname_chiral_field_final1,"fieldP_process_real.dat");}

        	chiral_field_final1= fopen(fname_chiral_field_final1,"r");
		
 
    		for (i=1; i<=NX; i++) {
      			for (j=1; j<=NY; j++) {
        			for (k=1; k<=NZ; k++) {
          						densityP[i][j][k] = 0.0;
          if (RERUN == 0){
	fscanf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",&temp,&temp,&temp,&psix[i][j][k],&psiy[i][j][k],&psiz[i][j][k]);
	//fscanf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",&temp,&temp,&temp,&psix[i][j][k],&psiy[i][j][k],&psiz[i][j][k],&modpsi[i][j][k]);
							fieldP[i][j][k] = 0.0;// This is field P - modpsi[i][j][k]/6.0;
		}
	else {
		fscanf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\n",&temp,&temp,&temp,&fieldP[i][j][k]); /*It is assumed that the fieldP has a value otained after subtracting the spatial average of the scalar field*/

		}
		
		}}}
			
			//if (RERUN == 0){				subtract_average(fieldP); /*Subtract spatial average of the scalar field*/
				//printf("%lf\t%lf\t%lf\t%lf\n",mx[i],my[j],mz[k],fieldP[i][j][k]);
		 	//		}
	
		 for (i=1; i<=NX; i++) {
                        for (j=1; j<=NY; j++) {
                                for (k=1; k<=NZ; k++) {
          						m = k + NZ*((j-1) + NY*(i-1));  
		      					/* For s=0 or time iteration = 0 and there are NITER + 1 iterations
		      					with index starting from 0 and reaching NITER on the last iteration*/

          						q_s[m] = 1.0;             /* Initial Condition for q(r,s)*/
          						qstar_s[m] = 1.0;         /* Initial Condition for qstar(r,s)*/

          						/* Spherical boundary condition */
          						if (mx[i]*mx[i]+my[j]*my[j]+mz[k]*mz[k] >= dimensions[3]*dimensions[3]) 
							{
            						q_s[m] = 0.0;
            						qstar_s[m] = 0.0;
          						}

          						x[m] = fieldP[i][j][k];
        						}
      						}
    					}
		fclose(chiral_field_final1);

    		/* Non-linear solver is called here */
    		ICOUNT = 0;
    		do {
      			iDONE = vecfun(g_total,d,x,out_f);
		  	printf("Step= %d\n",ICOUNT);
      
      		/*************** Write densities to file******************/
      		if (ICOUNT%100==0) {
        				FILE *densityP_process;
        				char fname_densityP_process[40];
        				sprintf(fname_densityP_process,"densityP_process_real.dat");
        				densityP_process = fopen(fname_densityP_process,"w");

        				fprintf(densityP_process,"%d\n",ICOUNT) ;
        				fprintf(densityP_process,"#x\t y\t z\t densityP\n") ;

					FILE *fieldP_process;
                                        char fname_fieldP_process[40];
                                        sprintf(fname_fieldP_process,"fieldP_process_real.dat");
                                        fieldP_process = fopen(fname_fieldP_process,"w");
                                        
        				for (i=1; i<=NX; i++) {
          					for (j=1; j<=NY; j++) {
            						for (k=1; k<=NZ; k++) {
              						fprintf(densityP_process,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k],densityP[i][j][k]) ;
              						fprintf(fieldP_process,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k],fieldP[i][j][k]) ;
            									}	
          								}
        							}
        				fclose(densityP_process);
        				fclose(fieldP_process);
      				} 
      		/****************/
  		free_energy = output_free();           
      		EvolveField(ICOUNT); //Evolution of field
		
      		ICOUNT = ICOUNT+1;

      		//iDONE = CheckFields();
    		}
    		while ((ICOUNT < LIMIT) && (!iDONE));//(ICOUNT<1);//
  
    		/********Output of the solver is distributed here ********/
    		for (i=1; i<=NX; i++) {
      			for (j=1; j<=NY; j++) {
        			for (k=1; k<=NZ; k++) {
          						m = k + NZ*((j-1) + NY*(i-1) );

          						fieldP[i][j][k] = x[m];
        						}
      						}
    					}
                   
    		/*******Write CONVERGED Monomer Density and Fields to a file ********/
    		FILE *monomer_density_final;
    		char fname_monomer_density_final[40];
    		sprintf(fname_monomer_density_final,"monomer_density_converged.dat");
                  
    		monomer_density_final = fopen(fname_monomer_density_final,"a");
               
    		FILE *neutral_field_final;
    		char fname_neutral_field_final[40];
    		sprintf(fname_neutral_field_final,"neutral_field_converged.dat");
    		neutral_field_final = fopen(fname_neutral_field_final,"a");

    		// fprintf(densityP_final,"%d\n",ICOUNT ) ;
    		for (i=1; i<=NX; i++) {
      			for (j=1; j<=NY; j++) {
        			for (k=1; k<=NZ; k++) {
          				fprintf(monomer_density_final,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k],densityP[i][j][k]) ; 
          				fprintf(neutral_field_final,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k],fieldP[i][j][k]) ; 
        						}
      						}
    					}
    		fclose(monomer_density_final);
    		fclose(neutral_field_final); 
                
 		return (free_energy);
}

int main() {
  		double free_energy_final;
  		int i,j,k;
  		clock_t begin,end;
  		double total_time;

  		begin = clock();

  		/*************Parameter Space**********************/
  		N = 100.0;    //strtod(argv[2],fp);   /* Degree of Polymerization*/
		kl = 1.0;                     /*Kuhn segment length */
  		RG = kl*sqrt(N/6.0);          /* Rg in nm*/
  		chi = 0.0;                    /* w_pp  */
  		n_anchor = 0; 	            /*No anchoring in this code*/
		para_D = chi/pow(RG,3.0);     /* D = chi/Rg^3  */
  		sigma = 100.0;         		/******Number of chains: NOT USED IN THIS CODE */
  		alpha_para = 4.4934;          /*Alpha parameter characterizing solution of the curl operator*/
  		Rsphere = 20.0;               /*Radius of sphere in units of RG*/
  		ds_step = 0.005;              /*****Chain contour step - ds*****/

		/*The following four parameters are used explicitly in this code*/ 
		Ke2 = 0.001;
		q0 = 0.0;
  		mu = 0.0;

  		NITER = (int)(1.0/ds_step);
  		//printf("%d\n",NITER);
  		initialize();   /* Allocation of memory*/
         
  		FILE *process_melt;
  		char fname_process_melt[40];

  		sprintf(fname_process_melt,"process_info_melt.dat");
  		process_melt = fopen(fname_process_melt,"a");
  		fprintf(process_melt,"N = %lf\t w = %lf\t wred = %lf\t sigma = %lf\n",N,chi,para_D,sigma);
  		fclose(process_melt);

  		/***********************************************/
  		/*****Dimensions of the cubic box (in units of RG)******/
  		dimensions[3] = Rsphere;        /*Radius of the sphere (in units of RG), which is equal to the 
                                    length of the cubic simulation box */
  		d[1] =  dimensions[3];
  		d[2] =  dimensions[3]; 
  		paramet[1] = para_D;
  		paramet[2] = sigma; 
  		free_energy_final = min_func_char(d,paramet);
        
  		/* Write the final results to the files*/
  		FILE *final_results;
  		char fname_final_results[40];
  		sprintf(fname_final_results,"final_results_melt.dat");
  		final_results = fopen(fname_final_results,"a");

  		fprintf(final_results,"N = %lf\t w = %lf\t wred = %lf\t sigma = %lf\n",N,chi,paramet[1],paramet[2]);
  		fprintf(final_results,"%e\t %e\t %e\t %e\n",free_energy_final,dimensions[1],dimensions[2],dimensions[3]) ;
  		fclose(final_results);
  
  		/*******Write CONVERGED Monomer Density and Fields to a file ********/
  		FILE *monomer_density_final1;
  		char fname_monomer_density_final1[40];
  		sprintf(fname_monomer_density_final1,"monomer_density_final.dat");
  		monomer_density_final1 = fopen(fname_monomer_density_final1,"a");

  		FILE *neutral_field_final1;
  		char fname_neutral_field_final1[40];
  		sprintf(fname_neutral_field_final1,"neutral_field_final.dat");
  		neutral_field_final1 = fopen(fname_neutral_field_final1,"a");

  		//FILE *chiral_field_final1;
  		//char fname_chiral_field_final1[40];
  		//sprintf(fname_chiral_field_final1,"chiral_field_final.dat");
  		//chiral_field_final1= fopen(fname_chiral_field_final1,"a");
  
  		//FILE *chiral_vector_final1;
  		//char fname_chiral_vector_final1[40];
  		//sprintf(fname_chiral_vector_final1,"chiral_vector_final.dat");
  		//chiral_vector_final1=fopen(fname_chiral_vector_final1,"a");

  		/*Calculate f vector**/
		calc_f_scft();  
  
		FILE *f_vector_final1;
  		char fname_f_vector_final1[40];
  		sprintf(fname_f_vector_final1,"f_vector_final.dat");
  		f_vector_final1=fopen(fname_f_vector_final1,"a");

  		for (i=1; i<=NX; i++) {
    			for (j=1; j<=NY; j++) {
      				for (k=1; k<=NZ; k++) {
        			// if(densityP[i][j][k]>0.00001){
        				fprintf(monomer_density_final1,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], densityP[i][j][k]);
        				fprintf(neutral_field_final1,"%lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], fieldP[i][j][k]);
        				//fprintf(chiral_field_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], psix[i][j][k],psiy[i][j][k],psiz[i][j][k]);
        				//fprintf(chiral_vector_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], tx[i][j][k],ty[i][j][k],tz[i][j][k]);	  
        				fprintf(f_vector_final1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",mx[i],my[j],mz[k], fx[i][j][k],fy[i][j][k],fz[i][j][k]);	  
      							}
    						}
    					fprintf(monomer_density_final1,"\n");
    					//fprintf(chiral_vector_final1,"\n");
    					fprintf(f_vector_final1,"\n");
  					}

  		fclose(monomer_density_final1);
  		fclose(neutral_field_final1);
  		//fclose(chiral_field_final1);
  		//fclose(chiral_vector_final1);
  		fclose(f_vector_final1);

  	

  		close_initialize();  
  
		end = clock();
  		total_time = (double)(end-begin);
  		printf("Total execution time= %lf\n",total_time);
  		return 0;
}
