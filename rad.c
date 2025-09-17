
enum{RHO,PPP,URR,UTH,UPH};

#include "header.h"

void open_h5( struct JetGrid * );
void close_h5( struct JetGrid * );
void freeGrid( struct JetGrid * );
void buildGrid( struct JetGrid * );
void printGrid( struct JetGrid * );
void fill_zone_data( int , int , int , double * , struct JetGrid * , int *);

int get_i( double , struct lightcurve * );
int get_wi( double , struct lightcurve * );
int get_wyi( double , struct lightcurve * );
double get_nu( int , struct lightcurve * );
double get_dtobs( int , struct lightcurve * );

double get_emissivity( double * theZone , struct lightcurve * lc , double nu , double t_c , struct par_list * theList ){

   double qem   = 1.75882009e+07;    // e/mc, 		Hertz/Gauss
   double qe3mc = 1.35351437e-22;    // e^3/mc^2,		ergs/Gauss
   double mp    = 1.67262178e-24;    // proton mass,	grams
   double me    = 9.10938291e-28;    // electron mass,	grams
   double sigma_T = 6.6524e-25;      // Thom. Xsection, cm^2

   double eps_B = theList->eps_B;
   double eps_E = theList->eps_E;
   double pp = theList->synch_p;
   double gam = theList->gamma_law;

   double rho0 = lc->rho0;
   double c = lc->c;
   //double tc = t_c*lc->r0/c;
   t_c *= lc->r0;
 
   double rho = theZone[RHO]*rho0;
   double Pp  = theZone[PPP]*rho0*c*c;

   //double X = theZone[4];
   double X = 1.0;
/*
   double wVol  = (gam-1.)*theZone[RHO] + theZone[PPP];
   double wCons = (gam-1.)*theZone[5]   + theZone[6];
   double epsRT = (wCons-wVol)/(theZone[PPP]);
   if( isnan(epsRT) || isinf(epsRT) ) epsRT = 0.0;
   if( epsRT < 1e-10 ) epsRT = 1e-10;
   if( epsRT > 1.0 ) epsRT = 1.0;
   eps_B = epsRT;
*/
   double n   = rho/mp;
   double e_th = Pp/(gam-1.);
   double B   = sqrt(8.*M_PI*eps_B*e_th);

   double gam_m = (pp-2.)/(pp-1.)*eps_E*(e_th)/(n*me*c*c);
   double gam_c = 6.*M_PI*me*c*c/( sigma_T * B*B * t_c );

   double nu_m  = (3./4./M_PI)*qem*gam_m*gam_m*B;
   double nu_c  = (3./4./M_PI)*qem*gam_c*gam_c*B;
 

   double Qnu = pow( nu/nu_m , 1./3. );
   if( nu > nu_m ) Qnu = pow( nu/nu_m , (1.-pp)/2. );

   if( theList->cooling ){
      double pow1 = 1./3.;
      double pow2 = (1.-pp)/2.;
      double pow3 = -pp/2.;
      double nu1 = nu_m;
      double nu2 = nu_c;
      //printf("nu_m/nu_c = %e\n",nu_m/nu_c); 
      if( nu_c < nu_m ){ nu1 = nu_c; nu2 = nu_m; pow2 = -.5; }

      if( nu < nu1 ){
         Qnu = pow( nu/nu1 , pow1 );
      }else if( nu < nu2 ){
         Qnu = pow( nu/nu1 , pow2 );
      }else{
         Qnu = pow( nu2/nu1 , pow2 )*pow( nu/nu2 , pow3 );
      }
   }

   double eps_m = (1./8./M_PI)*(pp-1.)*sqrt(3.)*qe3mc*n*B;
   //printf("Here in get emmisivity, eps_m = %e, Qnu = %e, return = %e\n",eps_m,Qnu,eps_m*Qnu);
   //printf("Here in get emmisivity, gam_m = %e, gam_c = %e, nu_m = %e, nu_c = %e\n",gam_m,gam_c,nu_m,nu_c);
   //printf("Here in get emmisivity, the zones are: theZone[0] = %e, theZone[1] = %e, theZone[2] = %e, theZone[3] = %e, theZone[4] = %e, theZone[5] = %e\n",theZone[0],theZone[1],theZone[2],theZone[3],theZone[4],theZone[5]);
   //printf("Here in get emmisivity, the zones are: theZone[0] = %e, theZone[1] = %e, theZone[RHO] = %e, theZone[PPP] = %e\n",theZone[0],theZone[1],theZone[RHO],theZone[PPP]);
   return( eps_m*Qnu*X );

}

double get_dV( double * theZone , int Nq , double r0 ){

   double rm  = theZone[Nq];
   double rp  = theZone[Nq+1];
   double thm = theZone[Nq+2];
   double thp = theZone[Nq+3];
   double phm = theZone[Nq+4];
   double php = theZone[Nq+5];

   double r = .5*(rp+rm);
   double dr = rp-rm;
   double th = .5*(thp+thm);
   double dth = thp-thm;
   double ph = .5*(php+phm);
   double dph = php-phm;

   double r2    = r*r;
   double sinth = sin(th);

   double dV = r2*sinth*dr*dth*dph;

   return( dV*r0*r0*r0 );

}

double get_doppler( double * theZone , int Nq , double th_obs ){

   double th = .5*(theZone[Nq+2] + theZone[Nq+3]);
   double ph = .5*(theZone[Nq+4] + theZone[Nq+5]);
 
   double ur  = theZone[URR];
   double uth = theZone[UTH];
   double uph = theZone[UPH];

   double gam = sqrt( 1. + ur*ur + uth*uth + uph*uph );
   double mu_r = cos(th)*cos(th_obs) + sin(th)*cos(ph)*sin(th_obs);
   double mu_th = cos(th)*cos(ph)*sin(th_obs) - sin(th)*cos(th_obs);
   //double vn = ( ur*cos(th)-uth*sin(th) )/gam;
   //double vn = sqrt(1.-(1./gam/gam));
   double vn = (ur*mu_r+uth*mu_th)/gam;
   //printf("phi = %e, theta = %e, ur = %e, uth = %e\n",ph,th,ur,uth);

   double dop = gam*(1.-vn);
   return( dop );

}

double time_dilate( double * theZone , double t ){
   double ur  = theZone[URR];
   double uth = theZone[UTH];
   double uph = theZone[UPH];
   double gam = sqrt( 1. + ur*ur + uth*uth + uph*uph);

   return( t/gam );
}

double get_tobs( double * theZone , int Nq , double t , double th_obs ){

   double rm  = theZone[Nq];
   double rp  = theZone[Nq+1];
   double thm = theZone[Nq+2];
   double thp = theZone[Nq+3];
   double phm = theZone[Nq+4];
   double php = theZone[Nq+5];
   double r = .5*(rm+rp);
   double th = .5*(thm+thp);
   double ph = .5*(phm+php);
   
   double mu_r = cos(th)*cos(th_obs) + sin(th)*cos(ph)*sin(th_obs);
   double mu_th = cos(th)*cos(ph)*sin(th_obs) - sin(th)*cos(th_obs);

   double t_obs = t - r*mu_r;
   //if(t_obs<1e-3 && t_obs>0.0){
	   //printf("t_obs = %e, t_obscal = %e, t = %e, r = %e, mu_r = %e\n", t_obs,t-r*mu_r,t,r,mu_r);}

   return( t_obs );

}

double get_SkyLoc( double * theZone , int Nq , double th_obs , double tobs ){

   double rm  = theZone[Nq];
   double rp  = theZone[Nq+1];
   double thm = theZone[Nq+2];
   double thp = theZone[Nq+3];
   double phm = theZone[Nq+4];
   double php = theZone[Nq+5];
   double r = .5*(rm+rp);
   double th = .5*(thm+thp);
   double ph = .5*(phm+php);
   
   double w1 = r*cos(th)*sin(th_obs);
   double w2 = r*sin(th)*cos(ph)*cos(th_obs);

   double w = w1 - w2;
   //printf("\nr = %e; th = %e; ph = %e; w = %e; t_obs = %e\n", r,th,ph,w,tobs);
   return( w );

}

double get_SkyLoc_Y( double * theZone , int Nq , double th_obs ){

   double rm  = theZone[Nq];
   double rp  = theZone[Nq+1];
   double thm = theZone[Nq+2];
   double thp = theZone[Nq+3];
   double phm = theZone[Nq+4];
   double php = theZone[Nq+5];
   double r = .5*(rm+rp);
   double th = .5*(thm+thp);
   double ph = .5*(phm+php);

   double wy = r*sin(th)*sin(ph);
   //printf("r = %e; th = %e; ph = %e; wy = %e\n", r, th, ph ,wy);

   return( wy );

}

void add_flux_lc( char * fname , struct lightcurve * lc , double dt , struct par_list * theList ){

   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   
   printf("From Rank = %d, Reading file %s...\n",rank,fname);
   struct JetGrid theGrid = {{0}};
   strcpy( theGrid.filename , fname );

   open_h5( &theGrid );
   buildGrid( &theGrid );
   printGrid( &theGrid );

   double t = theGrid.t;
   double dt_obs = 0.01;
   int Num_nu = lc->Nt;
   printf("Num_nu = %d\n",Num_nu);
   int Nq = theGrid.Nq;
   double theZone[Nq+6];	//Nq+rm,rp,tm,tp,pm,pp
   int i,j,k,jj,kk,frq_i;
   int Nt = theGrid.Nt;
   int Np = theGrid.Np;
   int * Nr = theGrid.Nr;
   int NtDiv = lc->Ndiv_th;
   int NpDiv = lc->Ndiv_ph;
   for( j=0 ; j<Nt ; ++j){
      for( k=0 ; k<Np ; ++k ){
         for( i=1 ; i<Nr[j] ; ++i ){
            for( jj=0 ; jj<NtDiv; ++jj ){
               for( kk=0 ; kk<NpDiv; ++kk){
                  int ZoneDiv[4] = {jj,kk,NtDiv,NpDiv};
                  fill_zone_data( i , j , k , theZone , &theGrid , ZoneDiv );
		  double rm  = theZone[Nq];
   		  double rp  = theZone[Nq+1];
   		  double thm = theZone[Nq+2];
   		  double thp = theZone[Nq+3];
   		  double phm = theZone[Nq+4];
   		  double php = theZone[Nq+5];
   		  double r = .5*(rm+rp);
   		  double th = .5*(thm+thp);
   		  double ph = .5*(phm+php);

                  double t_obs = get_tobs( theZone , Nq , t , lc->th_obs );
                  double w = get_SkyLoc(theZone , Nq , lc->th_obs, t_obs);  //get sky length w.r.t observer  
                  double wy = get_SkyLoc_Y(theZone , Nq , lc->th_obs);  //get sky height w.r.t observer
		  if(theList->to_do == 0){                  
                  int iobs = get_i( t_obs , lc );
                  //printf("t_obs = %e, i = %d\n", t_obs,iobs);
                  if( iobs > -1 && iobs < lc->Nt ){
                     double dt_obs = get_dtobs( iobs , lc );
                     double dop = get_doppler( theZone , Nq , lc->th_obs );
                     double t_c = time_dilate( theZone , t );
                     double eps = get_emissivity( theZone , lc , dop*(lc->nu) , t_c , theList );
                     double dV = get_dV( theZone , Nq , lc->r0 );
                     //printf("dV = %e eps = %e dt = %e dt_obs = %e\n",dV,eps,dt,dt_obs);
                     lc->F[iobs] += dV*eps*dt/dt_obs/dop/dop;
                        }
                      }
                  else if(theList->to_do == 1){
                     //printf("Am I here in to do????? Rank = %d, file = %s\n",rank,fname);
                     //printf("t_obs = %e, t_up = %e, t_lo = %e\n", t_obs,(lc->ti)+dt_obs,(lc->ti)-dt_obs);
                     if( (t_obs < (lc->ti)+dt_obs) && (t_obs > (lc->ti)-dt_obs) ){
                     //printf("TIME SATISFIED!!!!!! Rank = %d, file = %s\n",rank,fname);
                     for( frq_i=0 ; frq_i<lc->Nt ; ++frq_i){
                       //printf("From Rank = %d, Reading file %s...\n",rank,fname);
                       //printf("IS THE nu LOOP THE ISSUE????? Num_nu = %d\n", Num_nu);
                  	double nu = get_nu(frq_i,lc);
                  	//printf("nu = %e, From Rank = %d, Reading file %s...\n",nu,rank,fname);
      	          	//int iobs = get_i( t_obs , lc );
               	//printf("t_obs = %e, i = %d\n", t_obs,iobs);
                	//if( iobs > -1 && iobs < lc->Nt ){
                     	//double dt_obs = get_dtobs( iobs , lc );
                     	double dop = get_doppler( theZone , Nq , lc->th_obs );
                     	double t_c = time_dilate( theZone , t );
                     	double eps = get_emissivity( theZone , lc , dop*nu , t_c , theList );
                     	double dV = get_dV( theZone , Nq , lc->r0 );
                       //printf("dV = %e eps = %e dt = %e dt_obs = %e\n",dV,eps,dt,dt_obs);
                     	lc->F[frq_i] += dV*eps*dt/dt_obs/dop/dop;
                     	  }                  
                      }
                  }
                  else if(theList->to_do == 2){
                    if( (t_obs < (lc->ti)+dt_obs) && (t_obs > (lc->ti)-dt_obs) ){
                        //printf("\nAm I here in to do????? Rank = %d, file = %s, w = %e\n",rank,fname,w);
                        int iobs = get_wi( w , lc );
			int iobsy = get_wyi( wy , lc );
                        //printf("\n iobs = %d \n",iobs);
                        if( iobs > -1 && iobs < lc->Nt && iobsy > -1 && iobsy < lc-Nt){
                           double dop = get_doppler( theZone , Nq , lc->th_obs );
                     	    double t_c = time_dilate( theZone , t );
                     	    double eps = get_emissivity( theZone , lc , dop*(lc->nu) , t_c , theList );
                     	    double dV = get_dV( theZone , Nq , lc->r0 );
                     	    //printf("\ndV = %e eps = %e dt = %e dt_obs = %e\n",dV,eps,dt,dt_obs);
                     	    //lc->F_sky[ lc->Nt * iobsy + iobs] += ((rp - rm)/4.*M_PI)*eps*dt/dt_obs/dop/dop;
			    lc->F_sky[ lc->Nt * iobsy + iobs] += dV*eps*dt/dt_obs/dop/dop;
			    //printf("\nr = %e; w = %e; wy = %e; wi = %d; wyj = %d;  ti = %e;  t_obs = %e; t_code = %e; F = %e\n", r,w,wy,iobs,iobsy,lc->ti,t_obs,t,lc->F_sky[lc->Nt * iobsy + iobs]);
                         }
                       }                     
                      }
                   else{
                       printf("\n****** WRONG INPUT to_do ******");
                       	break; 
                    		} 
                    }
                  }
               }
            }
         }

   close_h5( &theGrid );
   freeGrid( &theGrid );   
}

double get_time( char * fname ){

   //printf("Getting time of %s...\n",fname);
   struct JetGrid theGrid = {{0}};
   strcpy( theGrid.filename , fname );

   open_h5( &theGrid );
   buildGrid( &theGrid );

   double t = theGrid.t;
   
   close_h5( &theGrid );
   freeGrid( &theGrid );

   //printf("t = %.2e\n",t);
   return(t);
}


