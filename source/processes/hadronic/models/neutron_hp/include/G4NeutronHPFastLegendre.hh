// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFastLegendre.hh,v 1.3 1999-07-02 09:58:57 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFastLegendre_h
#define G4NeutronHPFastLegendre_h 1

#include "globals.hh"

class G4NeutronHPFastLegendre
{
  public:
  
  G4NeutronHPFastLegendre()
  {
    value = new G4double * [31];
    value[0] = l0;
    value[1] = l1;
    value[2] = l2;
    value[3] = l3;
    value[4] = l4;
    value[5] = l5;
    value[6] = l6;
    value[7] = l7;
    value[8] = l8;
    value[9] = l9;
    value[10] = l10;
    value[11] = l11;
    value[12] = l12;
    value[13] = l13;
    value[14] = l14;
    value[15] = l15;
    value[16] = l16;
    value[17] = l17;
    value[18] = l18;
    value[19] = l19;
    value[20] = l20;
    value[21] = l21;
    value[22] = l22;
    value[23] = l23;
    value[24] = l24;
    value[25] = l25;
    value[26] = l26;
    value[27] = l27;
    value[28] = l28;
    value[29] = l29;
    value[30] = l30;
    integral = new G4double * [31];
    integral[0] = i0;
    integral[1] = i1;
    integral[2] = i2;
    integral[3] = i3;
    integral[4] = i4;
    integral[5] = i5;
    integral[6] = i6;
    integral[7] = i7;
    integral[8] = i8;
    integral[9] = i9;
    integral[10] = i10;
    integral[11] = i11;
    integral[12] = i12;
    integral[13] = i13;
    integral[14] = i14;
    integral[15] = i15;
    integral[16] = i16;
    integral[17] = i17;
    integral[18] = i18;
    integral[19] = i19;
    integral[20] = i20;
    integral[21] = i21;
    integral[22] = i22;
    integral[23] = i23;
    integral[24] = i24;
    integral[25] = i25;
    integral[26] = i26;
    integral[27] = i27;
    integral[28] = i28;
    integral[29] = i29;
    integral[30] = i30;
    
    G4int i;
    for(i=0;i<31;i++) theNbin[i]=1+200*(i+1);
  }
  
  ~G4NeutronHPFastLegendre()
  {
    delete [] value;
    delete [] integral;
  }
  
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
