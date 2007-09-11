

#ifndef CommonWrapper_hh
#define CommonWrapper_hh 1

// Constants:
#define FSIZEF 15
#define IGRAINESIZEF 19
#define MATSIZEF 500
#define MATGEOSIZEF 6
#define LGNSIZEF 9
#define LNSIZEF 30
#define SAXWROWSF 30 
#define SAXWCOLSF 500
#define DTONSIZEF 13
#define SPL2SIZEF 100
#define BL1SIZEF 300
#define BL2CROISSIZEF 19900
#define BL2INDSIZEF 19900
#define BL3SIZEF 300
#define BL5SIZEF 300
#define BL9SIZEF 300
#define VARSIZEF 3
#define VAEPSSIZEF 250
#define VAAVM 1000

extern "C" {
  struct hazard_common {
    int ial;	
    int IY[19];
  };
  extern struct hazard_common hazard_; 

  struct mat_common {
    int zmat[MATSIZEF];
    int amat[MATSIZEF];
    float bmax_geo[MATSIZEF][MATGEOSIZEF];
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
    float x[SAXWCOLSF][SAXWROWSF];
    float y[SAXWCOLSF][SAXWROWSF];
    float s[SAXWCOLSF][SAXWROWSF];
    int n;
    int k;
  };
  extern struct saxw_common saxw_;

  struct light_gaus_nuc_common {
    float rms1t[LGNSIZEF];
    float pf1t[LGNSIZEF];
    float pfln[LGNSIZEF];
    float tfln[LGNSIZEF];
    float vnuc[LGNSIZEF];
  };
  extern struct light_gaus_nuc_common light_gaus_nuc__;

  struct light_nuc_common {
    float r[LNSIZEF];
    float a[LNSIZEF];
  };
  extern struct light_nuc_common light_nuc__;

  struct calincl_common {
    float f[15];
    int icoup;
  };
  extern struct calincl_common calincl_;

//     float xsp[SPL2SIZEF];
//     float ysp[SPL2SIZEF];
//     float a[SPL2SIZEF];
//     float b[SPL2SIZEF];
//     float cc[SPL2SIZEF];
//     float nbp;

  struct spl2_common {
    float x[SPL2SIZEF];
    float y[SPL2SIZEF];
    float a[SPL2SIZEF];
    float b[SPL2SIZEF];
    float c[SPL2SIZEF];
    int n;
  };
  extern spl2_common spl2_;

  struct dton_common {
    double c[DTONSIZEF];
    double d[DTONSIZEF];
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

