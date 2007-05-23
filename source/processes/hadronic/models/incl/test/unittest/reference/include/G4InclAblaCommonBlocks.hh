
#ifndef InclAblaCommonBlocks_hh
#define InclAblaCommonBlocks_hh 1

// For debugging

extern "C" {
  extern struct {
    int cppevent;
    int allevents;
    int ntranscoul;
  } debugval_;
}

// C Dialogue with INCL for nucleus density and parameters.
//       COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX

extern "C" {
  extern struct {
    float R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX;
  } ws_;
}

// C Input for INCL, to be known at the CALL PNUC (SUPPRIMER icoup!!)
//       COMMON/CALINCL/FINPUT(15),icoup

extern "C" {
  extern struct {
    float FINPUT[15];
    int icoup;
  } calincl_;
}

// C ial generateur pour le cascade (et les IY pour eviter les correlations)
//       DIMENSION IY(19),IYV(19)
//       COMMON/hazard/ial,IY

extern "C" {
  extern struct {
    int ial;
    int IY[19];
  } hazard_;
}

// 	parameter (max=250)
// 	REAL*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS
// 	INTEGER AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS
//         COMMON/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
//      +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,FISPRO,NTRACK,
//      +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
//      +TETLAB(max),PHILAB(max)
#define VARNTP_MAX 250
extern "C" {
  extern struct {
    int MASSINI, MZINI;
    float EXINI;
    int MULNCASC, MULNEVAP,MULNTOT;
    float BIMPACT;
    int JREMN,KFIS;
    float ESTFIS;
    int IZFIS,IAFIS;
    float FISPRO;
    int NTRACK,ITYPCASC[VARNTP_MAX],AVV[VARNTP_MAX], ZVV[VARNTP_MAX];
    float  ENERJ[VARNTP_MAX],PLAB[VARNTP_MAX],TETLAB[VARNTP_MAX],PHILAB[VARNTP_MAX];
  } varntp_;
}

//       REAL*4 Bavat,TIME,ENERGY,EPSd,EPS2,EPS4,EPS6,EPSf
//       INTEGER Bloc_Paul,Bloc_CDPP,GO_OUT,avm,DEL1,DEL2
//       PARAMETER (avm=1000)
//       COMMON/VAR_AVAT/Kveux,Bavat,NOPART,NCOL,
//      s R1_in(3),R1_first_avat(3),
//      s EPSd(250),EPS2(250),EPS4(250),EPS6(250),EPSf(250),
//      s NB_AVAT,
//      s TIME(avm),L1(avm),L2(avm),JPARTL1(avm),JPARTL2(avm),
//      s DEL1(avm),DEL2(avm),ENERGY(avm),Bloc_Paul(avm),
//      s Bloc_CDPP(avm),GO_OUT(avm)
#define AVM 1000
extern "C" {
  extern struct {
    float Kveux,Bavat;
    int NOPART, NCOL;
    float R1_in[3],R1_first_avat[3];
    float EPSd[250],EPS2[250],EPS4[250],EPS6[250],EPSf[250];
    int NB_AVAT; //:::CHECK:::
    float TIME[AVM], L1[AVM], L2[AVM], JPARTL1[AVM], JPARTL2[AVM]; //:::CHECK:::
    int DEL1[AVM],DEL2[AVM];
    float ENERGY[AVM];
    int Bloc_Paul[AVM],Bloc_CDPP[AVM],GO_OUT[AVM];
  } varavat_;
}

extern "C" {
  extern struct {
    double probatrans;
    double randomtrans;
  } coulomb_;
}

#endif

