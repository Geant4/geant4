// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFastLegendre.hh,v 1.2 1999-06-29 18:43:58 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFastLegendre_h
#define G4NeutronHPFastLegendre_h 1

#include "globals.hh"

class G4NeutronHPFastLegendre
{
  public:
  
  G4NeutronHPFastLegendre();
  
  ~G4NeutronHPFastLegendre();
  
  G4double Integrate(G4int l, G4double costh)
  {
    G4int bin = GetBin(l, costh);
    G4double y1, y2;
//    G4cout <<"Testhpw G4NeutronHPFastLegendre::Integrate "<<l<<" "<<bin<<endl;
    y1 = integral[l][bin];
    y2 = integral[l][bin+1];
//    G4cout <<"Testhpw G4NeutronHPFastLegendre::Integrate exit"<<endl;
    return Interpolate(bin, l, y1, y2, costh);
  }
  
  inline G4double Evaluate(G4int l, G4double costh)
  {
    G4double result;
    G4int bin = GetBin(l, costh);
    if(bin != theNbin[l]-1)
    {
      G4double y1, y2;
      y1 = value[l][bin];
      y2 = value[l][bin+1];
      result = Interpolate(bin, l, y1, y2, costh);
    }
    else
    {
      result = value[l][bin];
    }
    return result;
  }
  private:
  
  inline G4int GetBin(G4int l, G4double costh)
  {
    G4int bin=0;
    bin = (theNbin[l]-1)*(costh+1)/2.;
    if(bin == theNbin[l]-1) bin--;
    return bin;
  }
  
  inline G4double Interpolate(G4int bin, G4int l, G4double y1, G4double y2, G4double x)
  {
    G4double slope = 0, off = 0, x1=0, x2=0, x1mx2;
    G4int half = (theNbin[l]-1)/2;
//    x1 = (bin-half)/G4double(half);
    x2 = (bin+1-half)/G4double(half);
    x1mx2 = 1./G4double( (theNbin[l]-1)/2 );
//    slope = (y2-y1)/(x2-x1);
    slope = (y2-y1)/x1mx2;
    off = y2-x2*slope;
    return x*slope+off;
  }
  
  G4double ** value;
  G4double ** integral;
  G4int theNbin[31];
  static G4double l0[201];
  static G4double i0[201];
  static G4double l1[401];
  static G4double i1[401];
  static G4double l2[601];
  static G4double i2[601];
  static G4double l3[801];
  static G4double i3[801];
  static G4double l4[1001];
  static G4double i4[1001];
  static G4double l5[1201];
  static G4double i5[1201];
  static G4double l6[1401];
  static G4double i6[1401];
  static G4double l7[1601];
  static G4double i7[1601];
  static G4double l8[1801];
  static G4double i8[1801];
  static G4double l9[2001];
  static G4double i9[2001];
  static G4double l10[2201];
  static G4double i10[2201];
  static G4double l11[2401];
  static G4double i11[2401];
  static G4double l12[2601];
  static G4double i12[2601];
  static G4double l13[2801];
  static G4double i13[2801];
  static G4double l14[3001];
  static G4double i14[3001];
  static G4double l15[3201];
  static G4double i15[3201];
  static G4double l16[3401];
  static G4double i16[3401];
  static G4double l17[3601];
  static G4double i17[3601];
  static G4double l18[3801];
  static G4double i18[3801];
  static G4double l19[4001];
  static G4double i19[4001];
  static G4double l20[4201];
  static G4double i20[4201];
  static G4double l21[4401];
  static G4double i21[4401];
  static G4double l22[4601];
  static G4double i22[4601];
  static G4double l23[4801];
  static G4double i23[4801];
  static G4double l24[5001];
  static G4double i24[5001];
  static G4double l25[5201];
  static G4double i25[5201];
  static G4double l26[5401];
  static G4double i26[5401];
  static G4double l27[5601];
  static G4double i27[5601];
  static G4double l28[5801];
  static G4double i28[5801];
  static G4double l29[6001];
  static G4double i29[6001];
  static G4double l30[6201];
  static G4double i30[6201];
};
#endif
