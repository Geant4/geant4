#ifndef G4BertiniCascade_h
#define G4BertiniCascade_h 1

#include "globals.hh"
#include "G4VIntraNuclearTransportModel.hh"

class G4BertiniModel;
class G4NucleusModel;
class G4BertiniData;

class G4BertiniCascade : public G4VIntraNuclearTransportModel {
public: // methods
  G4BertiniCascade();
  ~G4BertiniCascade();
  G4int operator == (G4BertiniCascade& right);
  G4int operator != (G4BertiniCascade& right);
  G4VParticleChange* ApplyYourself (const G4Track& aTrack, G4Nucleus& theNucleus);
  G4ReactionProductVector* Propagate (G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);
  void test();
private:
  void go();
  G4int chanwt();
  void relas(G4double minEnergy);
  void gene(G4double *z);
  void range();
  void recoil();
  void gthsig(G4int isgnl);
  void sgm(const G4int nsgnl, G4int it, const G4int is, const G4double de, 
const G4double em, G4int &i, G4double &e, G4double &s);
  void readh(G4double *geosig);
  void rainge(G4double eng, G4int mat, G4int ityp);
  void ck(G4int k, G4int npts, G4double *energy, G4double elas, G4double minEnergy);
  void hfield(G4double x, G4double y, G4double z, G4double u, G4double v, G4double w, G4double tip, G4double e);
  void getflt();
  void azirn(G4double &s, G4double &c);
  G4double adel(G4double e, G4double sum2, G4double sum1);
  G4double dxde(G4double ee, G4double sum1, G4double sum2);
  G4double exprnf();
  G4double gaurn();
  void getrig();
  G4double zfoi(const G4double z);
  G4double energyAZ(const G4double a, const G4double z);
  G4double enrg(const G4double a, const G4double z);
  void shxd(); // initialize data matrices locx and etc 
  G4double xlamb(const G4double x, const G4double y, const G4double z);
  void main();
  void main1();
  G4int main2();
  G4bool main3();
  void main4();
  void main5();
  void main6();
  void cascad();
  void update();
  void p1cli();
  G4double xsec(G4int id, G4int itp, G4double ec);
  void mainBodyInit();
  void mainBody();
  void mainBody2();
  void mainBody3();
  void mainBody3a();
  void mainBody3b();
  void mainBody3c();
  void mainBody4();
  void mainBody5();
  void mainBody6();
  void mainBody7();
  void mainBody8();
  void mainBody9();
  void mainFinal();
  void rou16(G4double *t);
  void input2();
  void datahi();
  G4double emin[6];
  void  scatt(G4int inc);
  void gtiso(G4double &u, G4double &v, G4double &w);	    
  G4int sprd(G4int ind, G4double rng, G4double rbtors, G4int mark, G4double dist);
  G4int gomprp(G4double x, G4double y, G4double z, G4double u, G4double v, G4double w, G4int nmed, 
G4double blz, G4double mark, G4double dist, G4double s, G4double xc, G4double yc, G4double zc, G4double , //:::
	       G4int mat, G4int nameol);
  /*
    void gomsor(G4int *x, 
    G4double *y, 
    G4double *z, 
    G4int nmed, 
    G4double blz );
  */
  G4double ecol(G4double rr, G4int mat, G4int ityp, G4double ec, G4double e);
  void modify(G4double &u, G4double &v, G4double &w, const G4double e, const G4double d,  const G4double tip);
  void dklos();
  void datalo();
  void erup();
  void hcol(G4int ib, G4int ityp, G4int hsig, G4int no1, G4double ec, G4int nopart, G4int *kind, G4double *ep,  
	    G4double *alf, G4double *bet, G4double *gam);
  void analz1(G4int i);
  void mfpd2(G4int i);
  void sors(G4int i);
  G4int iclock();
  G4double frnd();
private: 
  static const G4int    MAXPART = 59999;
  static const G4int    MAXP = 349;    
  enum MODEL {pseudo, elastic, inelastic, hydrogen, isobar, scaling};
  enum PARTICLE {eProton = 0, eNeutron = 1};
 
  G4NucleusModel *theTargetNucleus;
  G4BertiniModel *theBertiniModel;
  G4BertiniData  *theBertiniData;
  G4int ihec;
  //main
  G4double tmel1;
  G4int    jtimoh;
  //main1
  G4int nameol;
  G4double sigdk;
  G4double sigdki;
  G4double dkwto;
  G4double weig;
  G4double wto;
  G4double weight;
  G4double wtno;
  G4double wtnew;
  G4double aam[15];
  G4int ineg;
  G4int i310; // fluka interface
  G4int i300; // fluka interface
  // dklos
  G4int ikind;
  G4double tauinv[183];
  G4double tau[110];
  // relas
  G4int nelstp;
  G4int noelas;
  G4double com;
  G4int nledit;
  // range
  G4double eion[1][17];
  G4double sig2[20][5];
  //scatt
  G4double itfas[2]; 
  G4double efas[2];
  G4double alpfas[2];
  G4double betfas[2];
  G4double gamfas[2];
  // readh
  G4double crsc[2400];
  G4double geosig[240];
  G4int locf1[21];
  G4int locsig[21];
  G4double es1[500]; // original name was es
  G4double sige[500];
  // getflt
  G4double totels;
  G4double eps;
  G4double savels;
  G4double sigmx[7][17];
  // energy
  G4double waps[250][11];
  // enrg
  G4double cam2[130];
  G4double cam3[200];
  // modify
  G4double bfield[3];
  G4int locx[4][4];
  G4double etc[4][4];
  // xlamb
  G4int idgb;
  // main2
  G4int ipcl;
  G4int iqheh;
  G4int ncasca;
  G4int npsm;
  // main4
  G4int narriv;
  G4int ifircl;
  // main5
  G4int igo;
  G4double rbar;
  // main6
  G4double coupi;
  // mainFinal
  G4int ktim;
  G4int jtim;
  G4int itim;
  G4int jtime;
  G4int jtimo;
  G4double xtim;
  G4double xtime;
  G4double xtimo;
  G4double edum;
  G4double udum;
  G4double xdum;
  G4double ydum;
  G4double zdum;
  G4double wtdum;
  G4double wttot;
  G4double cneut;
  G4double vdum;
  G4double coupav;
  G4double nspltt;
  G4double cnhist;
  G4double cnhis;
  G4double edo5r;
  G4int icemwr;
  //xsec
  G4double hadronCrossSection[29850];
  // scatt
  G4int in;
  //  G4double ginum[183]; //::: solve duplicate name
  // hcol
  G4double e1nc;
  G4double dout[40]; //::: part of isobar common <- organize this common data together
  G4double alf[MAXP];
  G4double bet[MAXP];
  G4double itt[MAXP];
  G4double eff[MAXP];
  G4double alpp[MAXP];
  G4double bett[MAXP];
  G4double gamm[MAXP];
  G4double wtt[MAXP];
  G4int iqh;
  G4int iqk;
  G4double cdd;
  G4double eloss1;
  G4int icoun1; 
  G4int icoun2; 
  G4int icoumx; 
  G4int ke;
  G4int i1;
  G4double randi[4];
  G4double sqnm;
  G4int ipec[13];
  G4double esps[481];
  G4double rands[4];
  G4double clsm;
  G4double pxyz[4];
  G4double rlke;
  G4double value1;
  G4double pt[16];
  G4int i3;
  G4double snt;
  G4int ippp;
  G4int nno;
  G4double cst;
  G4double crdt[2];
  G4double rcpmv;
  // p1cli
  G4double p1oe1;
  // rou16
  G4int i2;
  // cascade
  G4double f[8];
  G4double prd[5]; 
  G4double pcap[10][17];
  G4double echek;
  G4double echekk;
  G4double echec;
  G4double npseub;  
  G4double npsub;
  G4double nozerb;
  G4double cobert;
  G4double cobcem;  
  G4double counsk;
  G4double coqevt;
  G4double npseus;
  G4double nozers;  
  G4int    nskol;
  G4double npseuq;
  G4double nskcol;
  G4int nok;
  G4double nozerq;
  G4double pdevst;  
  G4int kok;
  G4int iqev;
  G4int ibrt;
  G4int iheh;
  G4int iskl;
  G4int icr;
  G4int iphev;
  G4double eno;
  G4double tipno;
  G4double ecno;
  G4int nhist1;
  G4int nooo;
  G4int nonow;
  G4int io;
  G4int ou;
  G4int ihie;
  G4double tipp;
  G4double ectipp;
  G4int negzpr;
  G4double echicut;
  G4double ehicut;
  G4double r;
  G4int mxmat;
  G4int lmx;
  G4double pcaps;
  G4int nel[7]; // ::: COMON.F: nel(17), count.f: nel(20, 9, 10)
  G4double flag;
  G4double dkwt;
  G4int nf;
  G4double den[9][16];
  G4double zz[9][15];
  G4double r0;
  G4double r1;
  G4double rat;
  G4double delsig;
  G4double flt;
  G4int hsigg[5][16];
  G4double hsig;
  G4int nofsk;
  G4int no1;
  G4int ihit;
  G4int ihcol;
  G4int ibert;
  G4int ibb;
  G4int ibbarr[MAXP];
  G4int iibar[182]; // f77: IDMAX8=183
  G4int inc;
  G4double edpcol;
  G4double entrr;
  G4double weni;
  G4double anucco;
  G4double znucco;
  G4double copcol;
  G4double coupsd;
  G4double coupih;
  G4double einit;
  G4double esum;
  G4double heht;
  G4double hehel;
  G4double helnel;
  G4double anumt;
  G4double denn;
  G4double dif;
  G4double coelas;
  G4int ikin;
  G4double eheh[4];
  G4double edqheh;
  G4double couheh;
  G4double entr;
  G4double epk[MAXP];
  G4double eloss;
  G4double nnn;
  G4double sigg[9][16];
  G4double rcut;
  G4double rdif;
  G4double elas;
  G4int nl;
  G4int noel[16];
  G4double sgels[9];
  G4int idd;
  G4int id[9][16];
  G4double rath;
  G4double edsc;
  G4double csthcm;
  G4double ab3f1;
  G4double fone;
  G4double a[9][16];
  G4double cos1;
  G4double am2;
  G4double e2;
  G4double wlab;
  G4double p1;
  G4double bta;
  G4double wtlab;
  G4double zta;
  G4double gma;
  G4double wcm;
  G4double e4cm;
  G4double e3cm;
  G4double p3cm;
  G4double pz3cm;
  G4double pz3;
  G4double e3;
  G4double e4;
  G4double p3;
  G4double csth;
  G4double snth;
  G4double csph;
  G4double snph;
  G4double andit;
  G4double ctofen;
  G4double j;
  G4double denlm;
  G4double invalid;
  G4double sigglm;
  G4double lele;
  G4double matt;
  G4double anucc;
  G4double znucc;
  G4double ehinn;
  G4double ctofe;
  G4double eskal;
  G4double eskale[4];
  G4double ihecc;
  G4double counpi;
  G4double anuma;
  G4double anumas;
  G4double cobc;
  G4double cobb;
  G4int nbogus;
  G4int nhistm;
  G4int imash;
  G4int ib;
  G4double amm;
  G4double ami;
  G4double conver;
  G4double edbrt;
  G4double esum2;
  G4double pi0;
  G4double sume;
  G4double prono;
  G4double pipos;
  G4double pineg;
  G4double chgpis;
  G4double fpt;
  G4int ia;
  G4int iz;
  G4int ik;
  G4int nskwt;
  G4double edqev;
  G4double edqevv;
  G4double edskl;
  G4double edskll;
  G4double rmfas;
  G4double exfas;
  G4double refass;
  G4double coss;
  G4double npsus;
  G4double esumw;
  G4double refas;
  G4double zress;
  G4double aress;
  G4double sumsi;
  G4int icp;
  G4double zzz;
  G4double izz;
  G4double yzz;
  G4double xzz;
  G4double aaa;
  G4int iaa;
  G4double yaa;
  G4double xaa;
  G4double mtim;
  G4int ittim;
  G4double convrr;
  G4double amscom;
  G4double amasdf;
  G4double amasw;
  G4double wapsen;
  G4double amasd;
  G4double elab;
  G4double elosss;
  G4double elos;
  G4double adiff;
  G4double zdiff;
  G4int nn;
  //update
  G4int nomaxp;
  G4int i151[5999];
  G4int namea[MAXP];
  G4int kinda[MAXP];
  G4int ibarra[MAXP];
  G4double wta[MAXP];
  G4double tipa[MAXP];
  G4double ea[MAXP];
  G4double ua[MAXP];
  G4double va[MAXP];
  G4double wa[MAXP];
  G4int mxcess;
  G4int mobch;
  G4double wtaa;
  // rainge
  G4int iou;
  G4double einc;
  G4double fact;
  G4double xbar;
  G4double erng[899];
  G4double rnge[899][14];
  G4double sigr2[899][14];
  G4double sigr;
  G4double rngelo[6][14];
  G4double sigelo[6][14];
  G4double npsg[1][175];
  G4double pipsg[2][125];
  G4double emx[3];
  // gthsig
  G4double rms[6];
  G4double rmsinv[6];
  G4double elowd;
  G4double elow;
  G4int it;
  G4int is;
  G4int isa[3];
  G4double crossSection[29849];
  G4int il;
  G4double el;
  G4double sl;
  G4double eh;
  G4double sh;
  G4double hsigmx[6];
  // mainBody9
  G4double up;
  G4double ehin;
  G4int kindk[MAXP];
  G4int icon;
  // number of
  G4int nodat; 
  G4int noevt;
  G4int nomaxn;
  G4int nomaxh;
  G4int nofask;
  G4int nobch;
  G4int nomax;
  G4int nocas;
  G4int no;
  G4int noo;
  G4int no5rca;
  // n
  G4int n; 
  G4int ncol;
  G4int nbertp;
  G4int nkey;
  G4int nhcasc;
  G4int negex;
  G4int ncas;
  G4int nsav;
  G4int nfirst;
  G4int nhcas;
  G4int npsmd;
  G4int nil;
  G4int niil;
  G4int name[MAXPART];
  G4int namax;
  G4int nhis;
  G4int nhist;
  G4int nhistt;
  G4int npidk;
  G4int nhstp;
  G4int nquit;
  G4int nhkey;
  G4int nnow;
  G4int num;
  // abov, belo
  G4int nbelo;
  G4int nabov;
  G4int nbeloo;
  G4int nabovo;
  // alha beta gamma
  G4double alpha[MAXP];
  G4double beta[MAXP];
  G4double gam[MAXP];
  G4int nbeta;
  G4int ngam;
  // evap
  G4int nevaph;
  G4int nevapl;
  G4int nerup;
  G4int nerupl;
  // number of particles
  G4int nopart;
  G4int noparto;
  G4int nseudo;
  G4int nparto;
  G4int nspred;
  G4int npart[5];
  // neutron
  G4int neutp;
  G4int neutno;
  //  G4file   fNeutronp; //:::
  G4double eneut;
  // nc
  G4int ncc;
  G4int ncco;
  G4int ncmc;
  G4int ncmosc;
  G4int nmed[MAXPART];
  G4int ncoutp[8]; 
  // wt
  G4double wtav;
  G4double wtcut;
  G4double wtlow;
  G4double wtn;
  G4double wtfas[MAXP];
  G4double oldwt;
  G4double wt[MAXPART];
  G4int nwtc;
  // dc
  G4double adc;
  G4double bdc;
  G4double gdc;
  G4double adcc;
  G4double bdcc;
  G4double gdcc;
  G4double sumdc;
  G4double ang1[99][2][1];
  G4double sinth;
  G4double costh;
  G4double sinphi;
  G4double cosphi;
  G4double cneuav;
  G4double edtotn;
  G4double rng;
  G4double q;
  G4double d;
  G4double dpr;
  G4double rr;
  G4int kk;
  G4double rp;
  G4double bodge;
  G4double dd;
  G4int lelem;
  G4int ll;
  G4int m;
  G4int jk;
  G4double uu;
  G4double sig;
  G4int kindi[MAXP];
  G4double countu ; 
  // e
  G4double emax;
  G4double ehipi;
  G4double cascaehin;
  G4double elop;
  G4double ep[MAXP];
  G4double epart[99][1];
  G4double ec[MAXPART];
  G4double erec;
  G4double ex;
  G4double e1;
  G4double ec1;
  // energy, x, y, z
  G4double energy[MAXPART];
  G4double x[MAXPART];
  G4double y[MAXPART];
  G4double z[MAXPART];
  G4double es;
  G4double xs;
  G4double ys;
  G4double zs;
  G4double u[MAXPART];
  G4double v[MAXPART];
  G4double w[MAXPART];
  G4double xc[MAXPART];
  G4double yc[MAXPART];
  G4double zc[MAXPART];
  G4double blz[MAXPART];
  G4int mat;
  G4int mmf;
  G4int mark;
  G4int medeq[14];
  G4double dist;
  G4double bold;
  G4double bbtors;
  G4double rbtors;
  G4double s;
  G4double denh[16]; 
  G4double tip[MAXPART]; 
  G4double tipstr; 
  G4int lowaz;
  G4double apro;
  G4double apr;
  G4double zpro;
  G4double zpr;
  // i
  G4int iq;
  G4int ic;
  G4int ih;
  G4int ityp;
  G4int ind;
  G4int ibertp;
  G4int iqhist;
  G4int iltemp;
  G4int ibtemp;
  G4int icno;
  G4int icem;
  // maximum
  G4int maxcas;
  G4int maxbch;
  G4int maxcso;
  G4int maxbco;
  // stor
  G4double gmstor[24999];  
  G4int kstor;
  // evt
  G4double pevts[15];
  G4double revts[15];
  G4double phevts[15];
  G4double rhevts[15];
  G4double pdevts[15];
  G4double peevts[15];
  G4double reevts[15];
  G4double hpevt[15];
  G4double hrevt[15];
  G4double hpevth[15];
  G4double hrevth[15];
  G4double hpdevt[15];
  G4double hevsum;
  G4double peiv[15]; 
  G4double barr[MAXPART]; 
  G4int kind[349];
  G4int kindo[MAXPART];
  G4int ijevnt[29]; 
  G4double hepart[99][3]; 
  G4double sf;
  G4double poac[19];
  G4double ppac[19];

  G4double de;  // energy difference for tabulated cross section data
  static G4double ppnp[9][2]; 
  static G4double enrgy[9];
  static G4double ginum[30];
  static G4double fli[13];
  static G4double particleMass[7];
  static G4double particleCharge[7];
  static G4double massNucleon;
  static G4double massPionCharged;
  static G4double massPionZero;   
  static G4double oneThird;
  static G4double twoThirds;
  static G4double fourThirds;

  G4int ifirst;
  G4double massParticle[4];
  G4double pmdx[6426]; 
  G4double pmdd[6426];
  G4double eth[4][4];
  G4double eoprs;
  G4double coqheh;
  G4int inpt;
};

#endif





