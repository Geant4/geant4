#ifndef G4BertiniInelasticCollision_h
#define G4BertiniInelasticCollision_h 1

#include "globals.hh"
#include "G4BertiniModel.hh"

class G4BertiniInelasticCollision : public G4BertiniModel {
public:
  G4BertiniInelasticCollision();
  ~G4BertiniInelasticCollision();
private:
static const G4int MAXPART = 59999;
  void bert(G4int ibert, G4double *finput, G4int nopart, G4int *kind, G4double *eray, 
G4double *aray, G4double *bray, G4double *gray);
  void abran(G4int);
  void alp19();
  void alp28();
  void alpha();
  void angid(); 
  void bb();
  void bg6b();
  void bg6c(G4int int1);
  void bg6ca(G4int k, G4int l);
  void ccpes();
  G4int coll(G4int m);
  G4int collm(G4int m);
  void crjab(G4int k1, G4double pp);
  void dfmax();
  void ecpl();
  void frmic(G4double &gpart);
  void geo();
  void idk(); 
  void isom();
  void nn();
  void p1clc();
  void partin();
  void pfmax();
  void pstor();
  void punp(); 
  void signex();
  void spac32(G4int i);
  void spisom();
  void stor(G4double p0, G4int i, G4double e, G4double a, G4double b, G4double g, G4double w, 
	     G4int j, G4double erem, G4double ps, G4double wmass, G4int itype);     
  void store();
  void stph(G4double *w);
  void stpl(G4double *w);
  void stpr();
  void xyi(G4int ii, G4int jj, G4int kk);
  void undis();
  void zero();
  void spcn();
private: // data members
  //bg6c
  G4double abz;
  // weight struct
  G4double wtav;
  G4int nf;
  // tarset struct
  G4double anucc;
  G4double znucc;
  G4double signl;
  G4double anucco;
  G4double znucco;
  G4double ec1;
  // run struct
  G4int ke;
  // rn struct
  G4double count[9];
  // numpert struct
  G4int numm;
  G4int num;
  G4int numm2;
  G4int nhist;
  G4int nodat;
  G4int noevt;
  G4int numa;
  G4int numb;
  // munpu struct
  G4double pcc[11];
  G4double ppnb[4];
  G4int nip[12];
  // joint struct
  G4int ibertp;
  G4int nbertp;
  // iowert struct
  G4int io;
  // dens struct
  G4int lelem;
  G4int mat;
  G4double denlm;
  G4double sigglm;
  // heav1 struct
  G4double nhsum[3];
  G4double ehsum[3];
  G4double cehsum[3];
  // inout struct
  G4int iqiq;
  G4int ioo;
  // struct
  G4int nrt;
  G4int ln;
  G4double andit;
  G4double amasno;
  G4double zee;
  G4double einc;
  G4double ctofe;
  G4double prtin;
  G4double trsym;
  G4int randi[3];
  G4double sqnm;
  G4double rcpmv;
  G4int ipec[12];
  G4int isw[12];
  G4int iout[5];
  G4int nor;
  G4int inpt;
  G4int no;
  G4int nmas;
  G4int in;
  G4int iv;
  G4int ip; 
  G4int medium; 
  G4int inc;
  G4int ifca;
  G4int ifcc;
  G4int not;
  G4int it;
  G4int ifc;
  G4int knot;
  G4int i1;
  G4int i2;
  G4int i3;
  G4int i4;
  G4int i5; 
  G4int i6;
  G4int i7;
  G4int ik;
  G4int ka;
  G4int itote;
  G4int itoti;
  G4int itot2;
  G4int itot3;
  G4int nwds;
  G4double esps[480];
  G4double pt[47];
  G4double plvc[960];
  G4double pgvc[440];
  G4double pnbc[4];
  G4double space[177];
  G4double crsc[2399];
  G4double hvn[2];
  G4double hvp[2];
  G4double awd[2];
  G4double fvnp[2];
  G4double vnvp[2];
  G4double pmac[2];
  G4double ppan[2];
  G4double thpn[2];
  G4double ffptfn[2];
  G4double tffn[2];
  G4double tffp[2];
  G4double cfepn[5];
  G4double fmpn[5];
  G4double fmax[6];
  G4double s[36];
  G4double coordinate[2]; 
  G4double dcos[2];
  G4double curr[10];
  G4double d[5];
  G4double ce[20];
  G4double wkrpn[5];
  G4double pxyz[11];
  G4double c[2];
  G4double eco[1];
  G4double col[23];
  G4double cc[11];
  G4double pnidk[22];
  G4double out[39];
  G4double fcn;
  G4double fcp;
  G4double pgcnt;
  G4double pacnt;
  G4double pecnt;
  G4double value2;
  G4double pppda;
  G4double ppmda;
  G4double ppnda;
  G4double ppnna;
  G4double clcfe;
  G4double value1;
  G4int rands[4];
  G4double ans;
  G4double begru;
  G4double frand;
  G4int erand[3];
  G4double ex;
  G4double sign;
  G4double clsm;
  G4double efrp;
  G4double efrn;
  G4double strkp;
  G4double rlke;
  G4double p1oe1;  //:: abs removed 
  G4double polc;
  G4double pols;
  G4double sopc;
  G4double sops;
  G4double p2;
  G4double any;
  G4double sn;
   G4double energy[MAXPART];
  G4double absec;
  G4double com;
  G4double snt;
  G4double cst;
  G4double com2;
  G4double value3;
  G4double univ;
  G4double univer;
  G4double com1;
  G4double ftr;
  G4double com4;
  G4double com3;
  G4double a;
  G4double ctofen;
  G4double tapcrs[29849];
  G4double unive;
  G4int idum1[1];
  G4double adum1[5];
  G4double prtinn;
  G4double trsymm;
  G4double adum2[4];
  G4int idum2[59];
  G4double adum3[2113];
  G4double acrsc[2399];
  G4double adum4[119];
  G4double adum5[168];
  G4double adum6[2];
  G4double adum7[30];
  // novalu struct
  G4int noo;
  G4int nnow;
  G4double e1[11];
  // xyinc struct
  G4double x;
  G4double y;
  G4double d1;
  G4double ppec[125];
  G4double pmec[125];
  G4double ppscl[116];
  G4double pmscl[116];  
  G4double ppdc[6424];
  G4double dmin[100];
  G4double fmxsn[160];
  G4double fmxdn[129];
  G4double fsln[175]; 
  G4double pnec[125];
  G4double pmxc[125]; 
  G4double pnnec[125]; 
  G4double pnscl[116]; 
  G4double pnnsl[116]; 
  G4double ecn[175];
  G4double pec[175];
  G4double pspcl[157]; 
  G4double spcln[158];
  G4double pdpcl[129];
  G4double dpcln[129]; 
  // ppdcl(378) 
  G4int ir[3];
  G4double expa;
  // frmic
  G4double gpart;
};

#endif

