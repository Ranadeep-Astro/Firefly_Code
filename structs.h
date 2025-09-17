
struct par_list{
   double E_iso, n_ism;
   double eps_B, eps_E, synch_p;
   double gamma_law;
   double luminosity_distance, redshift, th_obs;
   char checkpoint_basename[256];
   char output_filename[256];
   int checkpoint_first;
   int checkpoint_last;
   int checkpoint_step;
   int Ndiv_theta, Ndiv_phi;
   double tobs_min;
   double tobs_max;
   double skySpan_min, skySpan_max;
   double skySpan_Y_min, skySpan_Y_max;
   int Num_Obs;
   double frequency_min;
   double frequency_max;
   int cooling;
   int to_do;
};

struct directory{

   int N;
   char ** fnames;
   double * t;
   double * dt;

};

struct lightcurve{

   double * F, * F_sky;
   double tm,tp,num,nup,skm,skp,skym,skyp;
   double nu,ti;
   int Nt, Ndiv_th,Ndiv_ph;
   double th_obs;
   double E0,rho0,r0,t0,c;
   double dL,z;

};

struct JetGrid{
   char filename[256];
   hid_t h5fil;
   hid_t grp1;
   hid_t grp2;
   double t;
   int Nt;
   int Np;
   int Nc;
   int Nq;
   int * Nr;
   int * Tindex;
   double * t_jph;
   double * p_kph;
};


