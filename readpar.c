enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "header.h"

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512];
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ) return(1);

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file( struct par_list * theList ){

   int rank,size;

   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

   char pfile[] = "in.par";

   int err=0;  

   int nrank;
   for( nrank=0 ; nrank<size ; ++nrank ){
      if( rank==nrank ){
         err += readvar( pfile , "E_iso"              , VAR_DOUB , &(theList->E_iso)              );
         err += readvar( pfile , "n_ism"              , VAR_DOUB , &(theList->n_ism)              );
         err += readvar( pfile , "eps_B"              , VAR_DOUB , &(theList->eps_B)              );
         err += readvar( pfile , "eps_E"              , VAR_DOUB , &(theList->eps_E)              );
         err += readvar( pfile , "synch_p"            , VAR_DOUB , &(theList->synch_p)            );
         err += readvar( pfile , "cooling"            , VAR_INT  , &(theList->cooling)            );
         err += readvar( pfile , "luminosity_distance", VAR_DOUB , &(theList->luminosity_distance));
         err += readvar( pfile , "redshift"           , VAR_DOUB , &(theList->redshift)           );
         err += readvar( pfile , "obs_theta"          , VAR_DOUB , &(theList->th_obs)             );
         err += readvar( pfile , "checkpoint_basename", VAR_STR  , &(theList->checkpoint_basename));
         err += readvar( pfile , "output_filename"    , VAR_STR  , &(theList->output_filename)    );
         err += readvar( pfile , "checkpoint_first"   , VAR_INT  , &(theList->checkpoint_first)   );
         err += readvar( pfile , "checkpoint_last"    , VAR_INT  , &(theList->checkpoint_last)    );
         err += readvar( pfile , "checkpoint_step"    , VAR_INT  , &(theList->checkpoint_step)    );
         err += readvar( pfile , "tobs_min"           , VAR_DOUB , &(theList->tobs_min)           );
         err += readvar( pfile , "tobs_max"           , VAR_DOUB , &(theList->tobs_max)           );
         err += readvar( pfile , "skySpan_min"        , VAR_DOUB , &(theList->skySpan_min)        );
         err += readvar( pfile , "skySpan_max"        , VAR_DOUB , &(theList->skySpan_max)        );
         err += readvar( pfile , "skySpan_Y_min"      , VAR_DOUB , &(theList->skySpan_Y_min)      );
         err += readvar( pfile , "skySpan_Y_max"      , VAR_DOUB , &(theList->skySpan_Y_max)      );
	 err += readvar( pfile , "Num_Obs"            , VAR_INT  , &(theList->Num_Obs)            );
         err += readvar( pfile , "frequency_min"      , VAR_DOUB , &(theList->frequency_min)      );
         err += readvar( pfile , "frequency_max"      , VAR_DOUB , &(theList->frequency_max)      );
         err += readvar( pfile , "Adiabatic_Index"    , VAR_DOUB , &(theList->gamma_law)          );
         err += readvar( pfile , "divide_theta"       , VAR_INT  , &(theList->Ndiv_theta)         );
         err += readvar( pfile , "divide_phi"         , VAR_INT  , &(theList->Ndiv_phi)           );
         err += readvar( pfile , "to_do"              , VAR_INT  , &(theList->to_do)              );
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   int errtot;
   MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );

   if( errtot > 0 ){
      printf("Read Failed\n");
      return(1);
   }

   return(0);

}


