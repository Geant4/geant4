#ifndef PYTHIA_WRAPPER_6_H
#define PYTHIA_WRAPPER_6_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, November 2000
// Version 6.200 update October 2001
// Wrapper for FORTRAN version of Pythia
// This wrapper is NOT intended as a part of HepMC - it is only supplied
// for your convenience.
//////////////////////////////////////////////////////////////////////////
//
// A simple example of calling Pythia from C++ using this header file is
// given in test/test_PythiaWrapper.cxx
//
// Note the pyhepc routine is used by Pythia to fill
// the HEPEVT common block uses double precision and 4000 entries.
//

#include <ctype.h>
#include <cstring>

//--------------------------------------------------------------------------
// Initialization routine


extern "C" {
  void initpydata(void);
}
#define initpydata initpydata_

//--------------------------------------------------------------------------
// PYTHIA Common Block Declarations

const int pyjets_maxn =4000;
extern "C" {
  extern struct {
    int n, npad, k[5][pyjets_maxn];
    double p[5][pyjets_maxn], v[5][pyjets_maxn];
  } pyjets_;
}
#define pyjets pyjets_

extern "C" {
  extern struct {
    int mstu[200];
    double paru[200];
    int mstj[200];
    double parj[200];
  } pydat1_;
}
#define pydat1 pydat1_

extern "C" {
  extern struct {
    int kchg[4][500];
    double pmas[4][500], parf[2000], vckm[4][4];
  } pydat2_;
}
#define pydat2 pydat2_

extern "C" {
  extern struct {
    int mdcy[3][500], mdme[2][8000];
    double brat[8000];
    int kfdp[5][8000];
  } pydat3_;
}
#define pydat3 pydat3_

extern "C" {
  extern struct {
    int mrpy[6];
    double rrpy[100];
  } pydatr_;
}
#define pydatr pydatr_

extern "C" {
  extern struct {
    int msel, mselpd, msub[500], kfin[81][2];
    double ckin[200];
  } pysubs_;
}
#define pysubs pysubs_

extern "C" {
  extern struct {
    int mstp[200];
    double parp[200];
    int msti[200];
    double pari[200];
  } pypars_;
}
#define pypars pypars_

extern "C" {
  extern struct {
    int mint[400];
    double vint[400];
  } pyint1_;
}
#define pyint1 pyint1_

extern "C" {
  extern struct {
    int iset[500], kfpr[2][500];
    double coef[20][500];
    int icol[2][4][40];
  } pyint2_;
}
#define pyint2 pyint2_

extern "C" {
  extern struct pin3 {
    double xsfx[81][2]; // Fortran is xsfx(2,-40:40)
    int isig[3][1000];
    double sigh[1000];
  } pyint3_;
}
#define pyint3 pyint3_

extern "C" {
  extern struct {
    int mwid[500];
    double wids[5][500];
  } pyint4_;
}
#define pyint4 pyint4_

extern "C" {
  extern struct pin5 {
    int ngenpd, ngen[3][501];   // Fortran is ngen(0:500,3)
    double xsec[3][501];        // Fortran is xsec(0:500,3)
  } pyint5_;
}
#define pyint5 pyint5_

extern "C" {
  extern struct pin7 {
    double sigt[6][7][7];   // Fortran is sigt(0:6,0:6,0:5)
  } pyint7_;
}
#define pyint7 pyint7_

extern "C" {
  extern struct pin8 {
    double xpvmd[13];   // Fortran is xpvmd(-6:6)
    double xpanl[13];   // Fortran is xpanl(-6:6)
    double xpanh[13];   // Fortran is xpanh(-6:6)
    double xpbeh[13];   // Fortran is xpbeh(-6:6)
    double xpdir[13];   // Fortran is xpdir(-6:6)
  } pyint8_;
}
#define pyint8 pyint8_

extern "C" {
  extern struct pin9 {
    double vxpvmd[13];  // Fortran is vxpvmd(-6:6)
    double vxpanl[13];  // Fortran is vxpanl(-6:6)
    double vxpanh[13];  // Fortran is vxpanh(-6:6)
    double vxpdgm[13];  // Fortran is vxpdgm(-6:6)
  } pyint9_;
}
#define pyint9 pyint9_

extern "C" {
  extern struct pssm {
    int imss[100];      // Fortran is imss(0:99)
    double rmss[100];   // Fortran is rmss(0:99)
  } pyssm_;
}
#define pyssm pyssm_

extern "C" {
  extern struct {
    double zmix[4][4];
    double umix[2][2];
    double vmix[2][2];
    double smz[4];
    double smw[2];
    double sfmix[4][16];
    double zmixi[4][4];
    double umixi[2][2];
    double vmixi[2][2];
  } pyssmt_;
}
#define pyssmt pyssmt_

extern "C" {
  extern struct {
    double rvlam[3][3][3];
    double rvlamp[3][3][3];
    double rvlamb[3][3][3];
  } pymsrv_;
}
#define pymsrv pymsrv_

extern "C" {
  extern struct prvnv {
    double ab[2][16][2];
    double rms[4];      // Fortran is rms(0:3)
    double res[5][6];
    int idr;
    int idr2;
    double dcmass;
    int kfr[3];
  } pyrvnv_;
}
#define pyrvnv pyrvnv_

extern "C" {
  extern struct prvpm {
    double rm[4];       // Fortran is rm(0:3)
    double a[2];
    double b[2];
    double resm[2];
    double resw[2];
    bool mflag;
  } pyrvpm_;
}
#define pyrvpm pyrvpm_

extern "C" {
  extern struct {
    double xxm[20];
  } pyints_;
}
#define pyints pyints_

extern "C" {
  extern struct {
    double x1;
  } pyg2dx_;
}
#define pyg2dx pyg2dx_

//--------------------------------------------------------------------------
// PYTHIA routines declaration

#define pyhepc pyhepc_
#define pyinit pyinit_
#define pylist pylist_
#define pystat pystat_
#define pyevnt pyevnt_
#define upinit upinit_
#define upevnt upevnt_
#define upveto upveto_
extern "C" {
  void pyhepc(int*);
  void pyinit(const char*,const char*,const char*,double*,int,int,int);
  void pylist(int*);
  void pystat(int*);
  void pyevnt();
  void upinit(){}
  void upevnt(){}
  void upveto(){}
}

// define methods to hide the subtle syntax necessary to call fortran from C++
inline void call_pyhepc( int mode ){ pyhepc( &mode ); }
inline void call_pyinit( const char* frame, const char* beam, const char* target,
                         double win )
{ pyinit( frame,beam,target,&win,strlen(frame),strlen(beam),strlen(target) ); }
inline void call_pylist( int mode ){ pylist( &mode ); }
inline void call_pystat( int mode ){ pystat( &mode ); }
inline void call_pyevnt(){ pyevnt(); }


//--------------------------------------------------------------------------
// PYTHIA block data
// ( with gcc it works to initialize the block data by calling
//   "pydata();" at beginning, but this fails for f77, so the fortran routine
//   initpydata.f is supplied ... call it instead for platform independent
//   behaviour )

#define pydata pydata_
extern "C" {
  void pydata(void);
}

#endif  // PYTHIA_WRAPPER_6_H
