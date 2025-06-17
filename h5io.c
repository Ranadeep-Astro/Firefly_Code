
#include "header.h"

void open_h5( struct JetGrid * theGrid ){
   hid_t h5fil = H5Fopen( theGrid->filename , H5F_ACC_RDWR , H5P_DEFAULT );
   theGrid->h5fil = h5fil;
   theGrid->grp1  = H5Gopen1( h5fil , "Grid");
   theGrid->grp2  = H5Gopen1( h5fil , "Data");
}

void close_h5( struct JetGrid * theGrid ){
   H5Gclose( theGrid->grp1  );
   H5Gclose( theGrid->grp2  );
   H5Fclose( theGrid->h5fil );
}

void getH5dims( hid_t h5grp , char * dset , hsize_t * dims ){
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
}

void readSimple( hid_t h5grp , char * dset , void * data , hid_t type ){
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
}

void readPatch( hid_t h5grp , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
}

void buildGrid( struct JetGrid * theGrid ){

   hsize_t dims[3];
   //Read in the time
   readSimple( theGrid->grp1 , (char *)"T" , &(theGrid->t) , H5T_NATIVE_DOUBLE );
   //Determine number of theta zones from the dimensions
   getH5dims( theGrid->grp1 , (char *)"t_jph" , dims );
   int Nt = dims[0]-1;
   getH5dims( theGrid->grp1 , (char *)"p_kph" , dims );
   int Np = dims[0]-1;

   //Allocate some crap based on that...
   theGrid->Nt = Nt;
   theGrid->Np = Np;
   theGrid->Nr = (int *) malloc( Nt*sizeof(int) );
   //printf("Nt = %d, Np = %d\n",Nt,Np);
   theGrid->Tindex = (int *) malloc( Nt*sizeof(int) );;
   theGrid->t_jph = (double *) malloc( (Nt+1)*sizeof(double) );
   theGrid->p_kph = (double *) malloc( (Np+1)*sizeof(double) );

   //Read in the Nt values of theta.
   readSimple( theGrid->grp1 , (char *)"t_jph" , theGrid->t_jph , H5T_NATIVE_DOUBLE );
   theGrid->t_jph++;
   //printf("t_jph = %e\n", * theGrid->t_jph);
   readSimple( theGrid->grp1 , (char *)"p_kph" , theGrid->p_kph , H5T_NATIVE_DOUBLE );
   theGrid->p_kph++;
   //printf("p_kph = %e\n", * theGrid->p_kph);

   int start[2]    = {0,0};
   int loc_size[2] = {Nt,Np};
   int glo_size[2] = {Nt,Np};

   //Read in structral data from the grid, to be used for the purposes of knowing where a given strip begins.
   readPatch( theGrid->grp1 , (char *)"Nr"    , theGrid->Nr     , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( theGrid->grp1 , (char *)"Index" , theGrid->Tindex , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);

   //Calculate the number of cells on the grid and the number of quantities per zone from the dimensions.
   getH5dims( theGrid->grp2 , (char *)"Cells" , dims );
   theGrid->Nc = dims[0];
   theGrid->Nq = dims[1]-1;

}

void freeGrid( struct JetGrid * theGrid ){
   if( theGrid->Nr    ) free( theGrid->Nr );
   if( theGrid->Tindex) free( theGrid->Tindex );
   if( theGrid->t_jph ){
      theGrid->t_jph--;
      free( theGrid->t_jph );
   }
}

void printGrid( struct JetGrid * theGrid ){
   printf("t = %.2e, Nt = %d, Np = %d\n",theGrid->t,theGrid->Nt,theGrid->Np);
   printf("Nc = %d Nq=%d\n",theGrid->Nc,theGrid->Nq);
}

void fill_zone_data( int i , int j, int k , double * theZone , struct JetGrid * theGrid , int * ZoneDiv ){

   int Nq = theGrid->Nq;
   int Nc = theGrid->Nc;

   int jDiv  = ZoneDiv[0];
   int kDiv  = ZoneDiv[1];
   int NjDiv = ZoneDiv[2];
   int NkDiv = ZoneDiv[3];

   double thDiv_m = (double)jDiv/(double)NjDiv;
   double thDiv_p = (double)(jDiv+1)/(double)NjDiv;
   double phDiv_m = (double)kDiv/(double)NkDiv;
   double phDiv_p = (double)(kDiv+1)/(double)NkDiv;

   int q;
   for( q=0 ; q<Nq+6 ; ++q ) theZone[q] = 0.0;

   int goodzone = (i>0 && j<theGrid->Nt);
   if( goodzone ) goodzone = i<theGrid->Nr[j];
   if( goodzone ){

      int start[2]    = { theGrid->Tindex[j] + i-1 , 0    };
      int loc_size[2] = { 2                        , Nq+1 };
      int glo_size[2] = { Nc                       , Nq+1 };
      double TrackData[ loc_size[0]*loc_size[1] ];
      readPatch( theGrid->grp2 , (char *)"Cells" , TrackData , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size);

      double rm = TrackData[Nq];
      double rp = TrackData[2*Nq+1];

      double t0 = theGrid->t_jph[j-1];
      double dt = theGrid->t_jph[j]-t0;

      double tm = t0 + thDiv_m*dt;//theGrid->t_jph[j-1];
      double tp = t0 + thDiv_p*dt;//theGrid->t_jph[j];
      
      double p0 = theGrid->p_kph[k-1];
      double dp = theGrid->p_kph[k]-p0;

      double pm = p0 + phDiv_m*dp;//theGrid->t_jph[j-1];
      double pp = p0 + phDiv_p*dp;//theGrid->t_jph[j];

      int q;
      for( q=0 ; q<Nq ; ++q ){
         theZone[q] = TrackData[Nq+1+q];
      }
      theZone[Nq]   = rm;
      theZone[Nq+1] = rp;
      theZone[Nq+2] = tm;
      theZone[Nq+3] = tp;
      theZone[Nq+4] = pm;
      theZone[Nq+5] = pp;

   }else{printf("Not good zone\n");}

}

void printZone( double * theZone , int Nq ){
   printf("rm = %e rp = %e\ntm = %e tp = %e\n",theZone[Nq],theZone[Nq+1],theZone[Nq+2],theZone[Nq+3]);
   printf("Primitive Variables = ");
   int q;
   for( q=0 ; q<Nq ; ++q ) printf(" %.2e",theZone[q]);
   printf("\n");
}



