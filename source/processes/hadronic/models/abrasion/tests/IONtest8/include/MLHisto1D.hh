#ifndef MLHisto1D_h
#define MLHisto1D_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "Histo1DVar.hh"
#include "CSVofstream.hh"
#include "RPTofstream.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLHisto1D : public Histo1DVar
{
public:
  MLHisto1D (G4String name,G4String unitx,G4String unity,
	     double *ep_list, size_t ep_list_len, side conv)
    :Histo1DVar (name, ep_list, ep_list_len, conv), theUnitx(unitx), theUnity(unity) {};
  MLHisto1D () : Histo1DVar(), theUnitx("[keV]"),theUnity("[counts/bin]")  {};
  ~MLHisto1D () {};

  void SetNormalisation (G4double dval) {nFac = dval;};
  G4double GetNormalisation () {return nFac;};
private:
  G4String theUnitx;
  G4String theUnity;
  G4double nFac;
  //
  //
  // INLINE DECLARATIONS/DEFINTIONS:
  //
public:
  inline void set_unitx (G4String new_unit) {theUnitx = new_unit;}
  inline void set_unity (G4String new_unit) {theUnity = new_unit;}
  //    inline const char *get_unit () {return theUnit.data();}
  inline const G4String get_unitx () {return theUnitx;}
  inline const G4String get_unity () {return theUnity;}
  //
  //
  // DECLARATIONS OF FRIELDS OF THE CLASS.
  //
  //  friend ostream & operator << (ostream &s, MLHisto1D &q);
  friend CSVofstream & operator << (CSVofstream &s, MLHisto1D &q);
  friend RPTofstream & operator << (RPTofstream &s, MLHisto1D &q);
};
////////////////////////////////////////////////////////////////////////////////
#endif
