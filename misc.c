
#include "header.h"

double get_time( char * );

void generate_directory( struct directory * dir , struct par_list * theList ){//, int rankd , int dsize ){

   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

   int Nmax  = theList->checkpoint_last;
   int Nmin  = theList->checkpoint_first;
   int Nskip = theList->checkpoint_step;

   int Nfiles = (Nmax-Nmin)/Nskip;
   int nmin = (Nfiles*rank)/size;
   int nmax = (Nfiles*(rank+1))/size;
   //printf("Size = %d; Rank = %d; Nfiles = %d; nmin = %d; nmax = %d\n", size,rank,Nfiles,nmin,nmax);

   if( rank >   0    ) --nmin;
   if( rank < size-1 ) ++nmax;
   //printf("Size = %d; Rank = %d; Nfiles = %d; nmin = %d; nmax = %d\n", size,rank,Nfiles,nmin,nmax);

   int N = nmax-nmin;
   dir->N = N;

   dir->fnames = (char **) malloc( N*sizeof( char * ) );
   int i;
   for( i=0 ; i<N ; ++i ) dir->fnames[i] = (char *) malloc( 256*sizeof(char) );
   dir->t  = (double *) malloc( N*sizeof(double) );
   dir->dt = (double *) malloc( N*sizeof(double) );

   int n;
   for( n = 0 ; n < N ; ++n ){
      int nn = Nmin+Nskip*(n+nmin);
      //printf("Rank = %d; nmin = %d; nmax = %d; n = %d; nn = %d\n", rank,nmin,nmax,n,nn);
      sprintf( dir->fnames[n] , "%s%04d.h5" , theList->checkpoint_basename , nn );
      if( n==N ) {
      MPI_Barrier( MPI_COMM_WORLD );
      printf("\nWhy am I here MIDWAY??\n");}
      dir->t[n] = get_time( dir->fnames[n] );
   }
   //MPI_Barrier( MPI_COMM_WORLD );
   //printf("\nI am here again End of Reading\n");
   double * t = dir->t;
   for( n = 0 ; n < N ; ++n ){
      double tm = t[n];
      if( n>0 ) tm = t[n-1];
      double tp = t[n];
      if( n < N-1 ) tp = t[n+1];
      double tc = t[n];
      double dtp = .5*(tp-tc);
      double dtm = .5*(tc-tm);
      dir->dt[n] = dtp+dtm;
   }

}

void free_directory( struct directory * dir ){
   int i;
   int N = dir->N;
   for( i=0 ; i<N ; ++i ) free( dir->fnames[i] );
   free( dir->fnames );
   free( dir->t  );
   free( dir->dt );
}

double get_t( int i , struct lightcurve * lc ){

   double Nt = (double)lc->Nt;
   double tm = lc->tm;
   double tp = lc->tp;

   double x = ((double)i+.5)/Nt;
   double t = pow( tp/tm , x ) * tm;

   return( t*(lc->t0) );

}

double get_w( int i , struct lightcurve * lc ){

   double Nt = (double)lc->Nt;
   double wm = lc->skm;
   double wp = lc->skp;

   double x = ((double)i+.5)/Nt;
   //double w = pow( wp/wm , x) * wm;
   double w = x*(wp - wm) + wm;

   return( w*(lc->r0) );

}


double get_wy( int i , struct lightcurve * lc ){

   double Nt = (double)lc->Nt;
   double wym = lc->skym;
   double wyp = lc->skyp;

   double x = ((double)i+.5)/Nt;
   //double wy = pow( wyp/wym , x) * wym;
   double wy = x*(wyp - wym) + wym;

   return( wy*(lc->r0) );

}

double get_nu( int i , struct lightcurve * lc ){

   //printf("Here in nu?????? Nt = %d\n", Num_nu);

   double Nt = (double)lc->Nt;
   //printf("WHAAAAT??????\n");
   double num = lc->num;
   double nup = lc->nup;
   //printf("Here in nu?????? A comeback!!!\n");

   double x = ((double)i+.5)/Nt;
   double nu = pow( nup/num , x ) * num;
   //printf("I don't think this is the issue, nu = %e Nt = %e\n",nu, Nt);

   return( nu );

}

double get_dtobs( int i , struct lightcurve * lc ){

   double Nt = (double)lc->Nt;
   double tm = lc->tm;
   double tp = lc->tp;

   double xm = ((double)i   )/Nt;
   double xp = ((double)i+1.)/Nt;
   double dt = pow( tp/tm , xp ) * tm - pow( tp/tm , xm ) * tm;

   return( dt );

}

int get_i( double t , struct lightcurve * lc ){

   double tm = lc->tm;
   double tp = lc->tp;
   double Nt = (double)lc->Nt;

   double x = log( t / tm ) / log( tp / tm );

   int i = (int) (Nt*x);
   if( t < 0.0 || x < 0.0 ) i = -1;
   if( x > 1.0 ) i = Nt;
   //printf("x = %e, i = %d\n", x,i);

   return(i);

}


int get_wi( double w , struct lightcurve * lc ){

   double wm = lc->skm;
   double wp = lc->skp;
   double Nt = (double)lc->Nt;

   double x = ((w-wm)/(wp-wm));
   //double x = log( w / wm) / log( wp / wm );

   int i = (int) (Nt*x);
   if( w < wm || x < 0.0 ) i = -1;
   if( w > wp || x > 1.0 ) i = Nt;
   //printf("\nHERE IN WI wm = %e, wp = %e, w = %e\n",wm,wp,w);
   //printf("x = %e, i = %d\n", x,i);

   return(i);

}

int get_wyi( double wy , struct lightcurve * lc ){

   double wym = lc->skym;
   double wyp = lc->skyp;
   double Nt = (double)lc->Nt;

   double x = ((wy-wym)/(wyp-wym));
   //double x = log( wy / wym) / log( wyp / wym );

   int i = (int) (Nt*x);
   if( wy < wym || x < 0.0 ) i = -1;
   if( wy > wyp || x > 1.0 ) i = Nt;
   //printf("\nHERE IN WI wm = %e, wp = %e, w = %e\n",wm,wp,w);
   //printf("x = %e, i = %d\n", x,i);

   return(i);

}

void initiate_lightcurve( struct lightcurve * lc , struct par_list * theList ){
//int Ntt , double tm , double tp , double nu , double E0 , double n0 , double dL , struct par_list * theList ){

   double Bethe = 1e51;			//erg
   double mp    = 1.67262178e-24;	//grams
   double c     = 2.99792458e10;	//cm/s
   double day   = 86400.;		//seconds
 
   double e0   = Bethe*theList->E_iso;
   double rho0 = mp*theList->n_ism;
   double r0 = pow( e0/rho0/c/c , 1./3. );
   double z = theList->redshift;
   double t0 = r0/c*(1.+z);
   char * fname = theList->output_filename;
   FILE * outFile = fopen(fname,"w");

   lc->c    = c;      				//cm/s
   lc->E0   = e0;      			//erg
   lc->rho0 = rho0;    			//g/cm^3
   lc->r0   = r0;      			//cm
   lc->t0   = t0/day;  			//days
   lc->z    = z;	      			//Gravitational Redshift

   int Nt = theList->Num_Obs;
   lc->Nt = Nt;
   lc->F = (double *) calloc( Nt , sizeof(double) );
   lc->F_sky = (double *) calloc( Nt * Nt , sizeof(double) );
   lc->tm = theList->tobs_min*day/t0;
   lc->tp = theList->tobs_max*day/t0;

   lc->num = theList->frequency_min*(1.+z);
   lc->nup = theList->frequency_max*(1.+z);
   lc->skm = theList->skySpan_min/r0;
   lc->skp = theList->skySpan_max/r0;
   lc->skym = theList->skySpan_Y_min/r0;
   lc->skyp = theList->skySpan_Y_max/r0;
   printf("\nskm = %e; skp = %e\n",lc->skm,lc->skp);
   lc->ti = theList->tobs_min*day/t0;
   lc->nu = theList->frequency_min*(1.+z);
   lc->dL = theList->luminosity_distance;
   lc->th_obs = theList->th_obs;
   if(theList->to_do == 0){
   printf("\n#########  CALCULATING LIGHTCURVE  #########\n");
   printf("\nFrequency = %e Hz\n", lc->nu);
   printf("\nStart time = %e sec\n",lc->tm*t0);
   printf("\nEnd time = %e sec\n",lc->tp*t0);
   fprintf(outFile,"#########  CALCULATING LIGHTCURVE  #########\n");
   fprintf(outFile,"Frequency = %e Hz\n", lc->nu);
   fprintf(outFile,"Start time = %e days\n",lc->tm);
   fprintf(outFile,"End time = %e days\n",lc->tp);
   fprintf(outFile,"E_iso = %e erg\n",lc->E0);
   fprintf(outFile,"n_ism = %e\n",theList->n_ism);
   fprintf(outFile,"Luminosity Distance = %e cm\n",lc->dL);
   fprintf(outFile,"Observer Angle = %e rad\n",lc->th_obs);
   fprintf(outFile,"epsilon_E = %e\n",theList->eps_E);
   fprintf(outFile,"epsilon_B = %e\n",theList->eps_B);
   fprintf(outFile,"p = %e\n",theList->synch_p);
   fprintf(outFile,"**********####################***********\n");
   fprintf(outFile,"\nTime(days)  Frequency(Hz)  Flux(mJy)\n");
   }
   else if(theList->to_do == 1){
   printf("\n#########  CALCULATING SPECTRUM  #########\n");
   printf("\nTime (seconds) = %e\n",lc->ti*t0);
   printf("\nStart Frequency (Hz) = %e\n", lc->num);
   printf("\nEnd Frequency (Hz) = %e\n", lc->nup);
   fprintf(outFile,"#########  CALCULATING SPECTRUM  #########\n");
   fprintf(outFile,"Time = %e days\n",lc->ti);
   fprintf(outFile,"Start Frequency = %e Hz\n", lc->num);
   fprintf(outFile,"End Frequency = %e Hz\n", lc->nup);
   fprintf(outFile,"E_iso = %e erg\n",lc->E0);
   fprintf(outFile,"n_ism = %e\n",theList->n_ism);
   fprintf(outFile,"Luminosity Distance = %e cm\n",lc->dL);
   fprintf(outFile,"Observer Angle = %e rad\n",lc->th_obs);
   fprintf(outFile,"epsilon_E = %e\n",theList->eps_E);
   fprintf(outFile,"epsilon_B = %e\n",theList->eps_B);
   fprintf(outFile,"p = %e\n",theList->synch_p);
   fprintf(outFile,"**********####################***********\n");
   fprintf(outFile,"\nTime(days)  Frequency(Hz)  Flux(mJy)\n");
   
   }
   else if(theList->to_do == 2){
   printf("\n#########  CALCULATING SKY MAP #########\n");
   printf("\nTime (seconds) = %e\n",lc->ti*t0);
   fprintf(outFile,"Frequency = %e Hz\n", lc->nu);
   fprintf(outFile,"#########  CALCULATING  MAP  #########\n");
   fprintf(outFile,"Time = %e days\n",lc->ti*lc->t0);
   fprintf(outFile,"Start Frequency = %e Hz\n", lc->num/(1.+z));
   fprintf(outFile,"End Frequency = %e Hz\n", lc->nup/(1.+z));
   fprintf(outFile,"E_iso = %e erg\n",lc->E0);
   fprintf(outFile,"n_ism = %e\n",theList->n_ism);
   fprintf(outFile,"Luminosity Distance = %e cm\n",lc->dL);
   fprintf(outFile,"Observer Angle = %e rad\n",lc->th_obs);
   fprintf(outFile,"epsilon_E = %e\n",theList->eps_E);
   fprintf(outFile,"epsilon_B = %e\n",theList->eps_B);
   fprintf(outFile,"p = %e\n",theList->synch_p);
   fprintf(outFile,"**********####################***********\n");
   fprintf(outFile,"\nTime(days)  Frequency(Hz)  Flux(mJy)  Sky Loc(cm)  Sky Height(cm)\n");
   
   }
   lc->Ndiv_th = theList->Ndiv_theta;
   lc->Ndiv_ph = theList->Ndiv_phi;
}

void free_lightcurve( struct lightcurve * lc ){
   if( lc->F ) free( lc->F );
   if( lc->F_sky ) free( lc-> F_sky );
}

void output_lc( char * fname , struct lightcurve * lc ){
   FILE * outFile = fopen(fname,"a");
   int Nt = lc->Nt;
   double z = lc->z;
   int i;
   double mJy = 1e26;
   double R2  = pow( lc->dL , 2. );
   for( i=0 ; i<Nt ; ++i ){
      double t = get_t( i , lc );
      fprintf(outFile,"%e %e %e\n",t,lc->nu,lc->F[i]*mJy/R2/(1.0+z));
   }
   fclose(outFile);
}

void output_spec( char * fname , struct lightcurve * lc ){
   FILE * outFile = fopen(fname,"a");
   int Nt = lc->Nt;
   double z = lc->z;
   int i;
   double t0 = lc->t0;			//Days
   double mJy = 1e26;			//miliJy
   double day   = 86400.;		//seconds
   double R2  = pow( lc->dL , 2. );
   for( i=0 ; i<Nt ; ++i ){
      double nu = get_nu( i, lc );
      fprintf(outFile,"%e %e %e\n",lc->ti*t0,nu,lc->F[i]*mJy/R2/(1.0+z));
   }
   fclose(outFile);
}

void output_skyLoc( char * fname , struct lightcurve * lc ){
   FILE * outFile = fopen(fname,"a");
   int Nt = lc->Nt;
   double z = lc->z;
   int i, j;
   double t0 = lc->t0;			//Days
   double r0 = lc->r0;			//cm
   double mJy = 1e26;			//miliJy
   double day   = 86400.;		//seconds
   double R2  = pow( lc->dL , 2. );
   //fprintf(outFile,"\n -----********************-----\n");
   for( i=0 ; i<Nt ; ++i ){
      double w = get_w( i, lc );      //cm
      for( j=0 ; j<Nt ; ++j){
        double wy = get_wy( j, lc );        //cm
        fprintf(outFile,"%e %e %e %e %e\n",lc->ti*t0,lc->nu,lc->F_sky[ Nt * j + i],w,wy);
	//int iobs = get_wi( w/r0 , lc );
        //int iobsy = get_wyi( wy/r0 , lc );
	//printf("\nWhile Printing...\nw = %e; wy = %e; wi = %d; i = %d; wyj = %d; j = %d;  ti = %e; F = %e\n", w,wy,iobs,i,iobsy,j,lc->ti,lc->F_sky[lc->Nt * j + i]);
      }
      }
   fclose(outFile);
}

