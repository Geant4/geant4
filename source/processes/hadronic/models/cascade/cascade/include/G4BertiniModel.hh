#ifndef G4BertiniModel_h
#define G4BertiniModel_h 1

#include "G4NucleusModel.hh"

#include "G4ReactionProductVector.hh"

#include "G4KineticTrackVector.hh"

class G4BertiniModel {
public:
  G4BertiniModel();
  ~G4BertiniModel();
  G4ReactionProductVector* Apply(G4KineticTrackVector* theSecondaries, G4NucleusModel* theTargetNucleus);
  void azio(G4double &s, G4double &c);
  G4double bovera(const G4double v, const G4double ve);
  void prot(G4double p0, G4double e, G4double p, G4double x, G4double y, G4double z);
  void pol1(G4double &x, G4double &y); 
  void signex(){};
  G4double exprn();
  void stor(G4double p0, G4int i, G4double e, G4double a, G4double b, G4double g,
	    G4double w, G4int j, G4double erem, G4double p, G4double wmass, G4int itype){};
  void spisom() {};
  void crdet(G4int nodata, G4double data[], G4double ener);
  void capunp() {};
  void castpr() {};
  void collm(G4int m) {};
  void cagene(G4double z) {};
  void coll(G4int m) {};
  void geo() {};
  void bg6ca(G4int, G4int) {};
  void crjab(G4int i, G4double t) {};
  void partin() {};
  void cacoll(G4int m) {};
  G4double dflran();
protected:
  static const G4int MAXPART = 59999; 
  static const G4int MAXP = 349;    
  // interpolateElasticNeutronData
  G4int nn;
  G4int idd;
  G4int locf1[21];
  G4int locsig[21];
  G4int noel[16];
  G4int id[9][16];
  G4double totels;
  G4double sgels[9];
  G4double ef[500];
  G4double f1[500];
  G4double es1[500]; //original name was es
  G4double sige[500];
  G4double den[9][16];
  G4double fone;
  // elastic collisions
  G4int ityp;
  G4double efas[MAXP];
  G4double alpfas[MAXP];
  G4double betfas[MAXP];
  G4double gamfas[MAXP]; 
  G4double wtfas[MAXP];
  G4double pb[4];
  G4double st0r;
  G4int itxxx[80]; // ::: fix name
  G4int numberOfNucleus;
  // rout1
  G4double out[39];
  G4double zee;
  G4double space[177];
  G4double amasno;
  G4double hvn[2];
  G4double hpw[2];
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
  G4double einc;
  G4double cfepn[5];
  G4double ctofen;
  G4double fmpn[5];
  G4double rands[3];
  G4double randi[3];
  G4double ppmda;
  G4double pppda;
  G4double ppnna;
  // rout2
  G4int i1;
  G4double value2;
  G4double fmax[6];
  G4double s[2];
  // rout3
  G4int isw[12];
  G4int no; 
  G4int inc;
  // rout4
  G4double curr[10];
  // rout5
  G4int not;
  // rout6
  G4double abz;
  G4int medium;
  G4int i3;
  G4double clsm;
  G4int knot;
  G4int it;
  G4double absec;
  // rout6a
  G4double strkp;
  G4int i2;
  G4double com;
  G4double rcpmv;
  G4double energy[MAXPART]; // ::: other e dimensions exists
  // rout7
  G4int ifca;
  G4double wkrpn[5];
  // rout7a
  G4double ex;
  // rout8
  G4int iv;
  G4int ifc;
  G4int ifcc;
  G4int in;
  G4double c[4];
  G4double d[6];
  // rou11
  G4int ik;
  G4double snt;
  // rou12
  G4double col[23];
  G4double ke;
  G4double pt[47];  
  // rou14
  G4double eco[1]; 
  G4double pnbc[4];
  G4double ccofe;
  G4double clcfe;
  G4double ctofe;
  // rou15
  G4double cst;
  // rou17
  G4double t[6426];
  G4double b[125];
  G4double rr[125]; //::: renamed from r
  G4double w[234]; //::: w[MAXPART] other sizes
  G4double g[233]; //::: hcol.f: g(55), bert1.f g(234)
  G4double rlke;
  G4double com2;
  G4double pacnt;
  // rou18
  G4double pnidk[22];
  G4double pxyz[3];
  // rout21
  G4double sqnm;
  G4double pgcnt;
  G4double univer;
  G4double pecnt;
  // rou22
  G4int i4;
  G4int i5;
  G4int i6;
  G4double coordinate[2];
  G4double dcos[2];
  G4double any;
  G4int ka;
  G4int nrt; 
  G4int ln;
  G4int ipec[12];
  G4int iout[5];
  G4int nor;
  G4int inpt;
  G4int nmas;
  G4int itote;
  G4int itoti;
  G4int itot2;
  G4int itot3;
  G4int nwds;
  G4int idum2[59];
  G4double esps[480];
  G4double plvc[960];
  G4double pgvc[440];
  G4double crsc[2399];
  G4double cc[11];
  G4double fcn;
  G4double fcp;
  G4double ppnda;
  G4double value1;
  G4int erand[3];
  G4int idum1[1];
  // numprt;
  G4int numm;
  G4int num;
  G4int numm2;
  G4int nhist;
  G4int nodat;
  G4int noevt;
  G4int numa;
  G4int numb;
  // nofaskk;
  G4int nofsk;
  // inpu
  G4double andt;
  G4double ctof;
  G4double ctofn;
  // joint 
  G4int ibertp;
  G4int nbertp;
  // isobar //::: eliminate ?
  G4int iswd1[15];
  G4int nno;
  G4int iiv;
  G4int ippp;
  G4int iswd2[4];
  //isob
  G4double nter;
  G4double cdd;
  G4double calcin[2];
  G4double value3;
  G4double tesiso;
  // count2
  G4int nolamb;
  G4int elamb;
  G4int nosigm;
  G4int esigm;
  G4int noprot;
  G4int eprot;
  G4int noneut;
  G4int eneutr;
  G4int nopion;
  G4int epion;
  G4int nopi0;
  G4int epi0;
  G4int nohigh;
  G4int ehigh;
  // xpd
  G4int npsg[1][175];	
  G4double pipsg[1][125];
  G4double hsigmx[6];
  G4int locx[3][3];
  G4double eth[3][3];
  // inout
  G4int intap;
  G4int ioutp;
  // ret
  G4int iretr;
  G4int idumno;
  G4int icoun1;
  G4int icoun2;
  G4int icoumx;
  G4int nopart;
  // gtdata
  G4double ginum[182];
  G4double gsqmas[182];
  G4int led[998];
  G4int ld[MAXP];
  G4int nct;
  G4int ncp;
  G4int ncl;
  // comon3
  G4int kindo[MAXPART];
  G4int kindi[MAXP];
  G4int kinda[MAXP];
  G4int kindb[MAXP];
  G4int ibbarr[MAXP];
  // weno
  G4int ibarra[MAXP];
  G4int ibarrb[MAXP];
  // iqhh
  G4int iqh;
  G4int iqk;
  // dcintp
  G4double z[175]; // :::ber1.f: z(126)
  G4int kind[MAXP];
  //xsec
  G4double sign;
  G4double univ;
  G4double unive;
  //qlp
  G4double r; //::: check rr / r
  G4double ec[MAXPART];
  G4int noo;
  G4int icon;
  G4double e1;
  G4double itfas[2]; 
  G4int verboseLevel;
  G4double massParticle[4];
  static G4double massNucleon;
  G4double casesn;
  G4double sf;
  G4double poac[19];
  G4double ppac[19];
  static G4double massPionCharged;
  static G4double massPionZero;
  static G4double oneThird;
  static G4double twoThirds;
  static G4double fourThirds;
  static G4double bindingEnergy;
  G4double crdt[24];     
  G4double begru;
  G4double p2;
  G4double ip;
  G4double pspcl[158];
  G4double pec[176];
  G4double spcln[158];
  G4double ecn[176];
  G4double cratio;
  G4double ppscl[117];
  G4double ppec[126];
  G4double pmscl[117];
  G4double pmec[126];
  G4double pmxc[126];
  G4double pnscl[117];
  G4double pnec[126];
  G4double pnnsl[117];
  G4double pnnec[126];
  G4double pdpcl[130];
  G4double dpcln[130];
  G4double fripn[117];
  G4double pnmi[101];
  G4double fmxsp[117];
  G4double pcfsl[234];
  G4double pnfsl[234];
  G4double dmin[101];
  G4double fmxsn[161];
  G4double fmxdn[130];
  G4double fsln[176];
  G4double sopc;
  G4double sops ;
  G4double frinn[161];
};

#endif






