#ifndef G4BertiniIsobarModel_h
#define G4BertiniIsobarModel_h 1

#include "G4BertiniModel.hh"

class G4BertiniIsobarModel : public G4BertiniModel
{
public:
  G4BertiniIsobarModel();
  ~G4BertiniIsobarModel();
private:
  void pol1(G4double& cs, G4double& si);
  void nucnuc(G4double p0, G4int nofas, G4int itype);
  void pinuc(G4double p0, G4int nofas, G4int itype);
  void prob(G4int itype, G4double p0);
  void rout1();
  void rout2(G4double *t);
  void rout3();
  void rout4();
  void rout5(G4double *t, G4double *b, G4double *rr);
  void rout6();
  void rout6a();
  void rout7();
  void rout7a();
  void rout8();
  void rou10();
  void rou11();
  void rou12();
  void rou13();
  void rou14();
  void rou15();
  void rou17();
  void rou18();
  void rou19();
  void rou20(G4double *dcin, G4double *dcln, G4double *dchn, G4double *pdci, G4double *pdch);
  void rou21(G4double *v, G4double *w, G4double *x, G4double *y, G4double *z);
  void rou22(G4double *v, G4double *w, G4double *x, G4double *y, G4double *z);
  void pi(G4double po, G4int itip, G4double ener, G4double pmom, G4double a, G4double b, G4double g);
  void q();
  void signex();
  void cole4(const G4double abz);
  G4double big7(G4double c, G4int ix);
  void dcpr(G4int lk);
  void dcintp(G4double *w);
  void pinst();
  G4double mud();
  void isob();
  void qrdet(G4int nodata, G4double *data, G4double ener);
  void qou17(G4double *fripn, G4double *pnmi, G4double *fmxsp, G4double *pcfsl, G4double *pnfsl);
  void qollm();
  void qou18(); 
  void qou21(G4double *frinn, G4double *dmin, G4double *fmxsn, G4double *fmxdn, G4double *fsln);   
  void qou19(); 
  void qlp19();
  void qlp28();
  void qlpha(); 
  void qngid(); 
  void qoll(); 
  void qene(G4double *z);
  void qdk();
  void qstor(); 
private: // data members
  static G4double angle[15]; // differential cross-section data 
  static G4double pmxxx[12]; // ::: pm renamed to pmxxx because of name clash
  static G4double dndpim[12][3]; // dn/dpi-
  static G4double dndpip[12][3]; // dn/dpi+
  static G4double cpimk[12][3];
  static G4double cpimnu[12][3];
  static G4double cpipk[12][3]; 
  static G4double cpipnu[12][3];
  static G4double cpk[12][3];
  static G4double cpnu[12][3];
  static G4double pz[12];
};

#endif



