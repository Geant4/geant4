
extern "C" {
  void ribm_(float *rndm, int *ial);

  void rgauss_(float *xg);

  void init_incl__(int *init_graine);
  void init_mat__(int *izmat, int *iamat, int *imat);
  void integ_(float *ami, float *ama, float *dr, float *fonc, float *res);
  void dens_deut__();
  void spl2ab_();

  float wsax_(float *r);
  float derivwsax_(float *r);
  float dmho_(float *r);
  float derivmho_(float *r);
  float derivgaus_(float *r);

  float flin_(float *xv); 
  void flin2_(int *k);
  float deutv_(int *l, float *q);
  float fm2_(int *j);
  float dens_(float *q);
  float splineab_(float *xv);

  double texp_(double *x);
  void sig_reac__(int *iprojo, double *E, double *A, float *sig);
  double radi_us__(double *a);

  void force_abs__(int *iprojo, float *at, float *zt, float *ep,float *bmax, float *pt, float *proba);
  void xabs2_(double *zp, double *ap, double *zt, double *at, double *ep, float *sig);
  void coulomb_transm__(float *E, float *fm1, float *z1, float *fm2, float *z2, float *proba);
  double clmb1_(double *ro, double *eta, double *t1);
  double clmb2_(double *ro, double *eta, double *t1);
  
  // Functions from CERNLIB
  float ranf_();
};
