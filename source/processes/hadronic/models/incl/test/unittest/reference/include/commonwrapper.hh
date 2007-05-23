

#ifndef CommonWrapper_hh
#define CommonWrapper_hh 1

// Constants:
#define FSIZE 15
#define IGRAINESIZE 19
#define MATSIZE 500
#define MATGEOSIZE 6
#define LGNSIZE 9
#define LNSIZE 30
#define SAXWROWS 30 
#define SAXWCOLS 500
#define DTONSIZE 13
#define SPL2SIZE 100
#define BL1SIZE 300
#define BL2CROISSIZE 19900
#define BL2INDSIZE 19900
#define BL3SIZE 300
#define BL5SIZE 300
#define BL9SIZE 300
#define VARSIZE 3
#define VAEPSSIZE 250
#define VAAVM 1000

extern "C" {
  struct hazard_common {
    int ial;	
    int IY[19];
  };
  extern struct hazard_common hazard_; 

  struct mat_common {
    int zmat[MATSIZE];
    int amat[MATSIZE];
    float bmax_geo[MATSIZE][MATGEOSIZE];
    int nbmat;
  };
  extern struct mat_common mat_;

  struct ws_common {
    float r0, adif, rmaxws, drws; 
    int nosurf; 
    float xfoisa; 
    int npaulstr; 
    float bmax;
  };
  extern struct ws_common ws_;

  struct saxw_common {
    float x[SAXWCOLS][SAXWROWS];
    float y[SAXWCOLS][SAXWROWS];
    float s[SAXWCOLS][SAXWROWS];
    int n;
    int k;
  };
  extern struct saxw_common saxw_;

  struct light_gaus_nuc_common {
    float rms1t[LGNSIZE];
    float pf1t[LGNSIZE];
    float pfln[LGNSIZE];
    float tfln[LGNSIZE];
    float vnuc[LGNSIZE];
  };
  extern struct light_gaus_nuc_common light_gaus_nuc__;

  struct light_nuc_common {
    float r[LNSIZE];
    float a[LNSIZE];
  };
  extern struct light_nuc_common light_nuc__;

  struct calincl_common {
    float f[15];
    int icoup;
  };
  extern struct calincl_common calincl_;

//     float xsp[SPL2SIZE];
//     float ysp[SPL2SIZE];
//     float a[SPL2SIZE];
//     float b[SPL2SIZE];
//     float cc[SPL2SIZE];
//     float nbp;

  struct spl2_common {
    float x[SPL2SIZE];
    float y[SPL2SIZE];
    float a[SPL2SIZE];
    float b[SPL2SIZE];
    float c[SPL2SIZE];
    int n;
  };
  extern spl2_common spl2_;

  struct dton_common {
    double c[DTONSIZE];
    double d[DTONSIZE];
    double fn;
  };
  extern dton_common dton_;
}

extern hazard_common *gHazard;
extern mat_common *gMat;
extern ws_common *gWs;
extern saxw_common *gSaxw;
extern light_gaus_nuc_common *gLight_gaus_nuc;
extern light_nuc_common *gLight_nuc;
extern calincl_common *gCalincl;
extern spl2_common *gSpl2;
extern dton_common *gDton;

#endif

