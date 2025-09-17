
#include "header.h"

void initiate_lightcurve( struct lightcurve * , struct par_list * );
void free_lightcurve( struct lightcurve * );

int read_par_file( struct par_list * );

void generate_directory( struct directory * , struct par_list * );
void free_directory( struct directory * );

void output_lc( char * , struct lightcurve * );
void output_spec( char * , struct lightcurve * );
void output_skyLoc( char * , struct lightcurve * );

void add_flux_lc( char * , struct lightcurve * , double , struct par_list * );
void add_flux_spec( char * , struct lightcurve * , double , struct par_list * );

int main( int argc , char * argv[] ){

   MPI_Init( &argc , &argv );
   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

   struct par_list theList;
   int err = read_par_file( &theList );
   if( err ) return(1);
   
   struct lightcurve lc = {0};
   initiate_lightcurve( &lc , &theList );
   int Nt = lc.Nt;
   
   /*else if(theList.to_do == 1){
   struct lightcurve lc = {0};
   initiate_lightcurve( &lc , &theList );
   int Nt = lc.Nt;}*/

   struct directory dir = {0};
   generate_directory( &dir , &theList );
   for(int n=0 ; n<dir.N ; ++n ){
   printf("Files loaded = %s\n",dir.fnames[n]);
   }
   MPI_Barrier( MPI_COMM_WORLD );
   //printf("\n!!!!Letting everyone on same page!!!!!\n");

   int n;
   for( n=0 ; n<dir.N ; ++n ){
      if( (n>0 || rank==0) && ( n<dir.N-1 || rank==size-1) )
         add_flux_lc( dir.fnames[n] , &lc , dir.dt[n] , &theList );
   	}

   double F[Nt];
   MPI_Allreduce( lc.F , F , Nt , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   //int i;
   //for( i=0 ; i<Nt ; ++i ) lc.F[i] = F[i];
   double F_sky[Nt * Nt];
   MPI_Allreduce( lc.F_sky , F_sky , Nt*Nt , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   int i,j;
   for( i=0 ; i<Nt ; ++i ){
      lc.F[i] = F[i];
      for( j=0 ; j<Nt ; ++j){
         lc.F_sky[ Nt * j + i] = F_sky[ Nt * j + i];
	   }
   }
 
   if( rank==0 && theList.to_do==0) output_lc( theList.output_filename , &lc );
   if( rank==0 && theList.to_do==1) output_spec( theList.output_filename , &lc );
   if( rank==0 && theList.to_do==2) output_skyLoc( theList.output_filename , &lc );

   free_lightcurve( &lc );
   free_directory( &dir );

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return(0);

}

