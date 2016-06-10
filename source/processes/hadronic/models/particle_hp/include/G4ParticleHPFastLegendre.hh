//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPFastLegendre_h
#define G4ParticleHPFastLegendre_h 1

#include "globals.hh"

class G4ParticleHPFastLegendre
{
  public:
  
  G4ParticleHPFastLegendre()
  {
    value = new const G4double * [31];
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
    integral = new const G4double * [31];
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
  
  ~G4ParticleHPFastLegendre()
  {
    delete [] value;
    delete [] integral;
  }
  
  G4double Integrate(G4int l, G4double costh)
  {
    if(l>30) return regularIntegrate(l,costh);
    G4int bin = GetBin(l, costh);
    G4double y1, y2;
//    G4cout <<"Testhpw G4ParticleHPFastLegendre::Integrate "<<l<<" "<<bin<<G4endl;
    y1 = integral[l][bin];
    y2 = integral[l][bin+1];
//    G4cout <<"Testhpw G4ParticleHPFastLegendre::Integrate exit"<<G4endl;
    return Interpolate(bin, l, y1, y2, costh);
  }
  
  inline G4double Evaluate(G4int l, G4double costh)
  {
    if(l>30) return regularEvaluate(l,costh);
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

      G4double regularEvaluate( int l , double x );
      G4double regularIntegrate( int l , double x );

  
  inline G4int GetBin(G4int l, G4double costh)
  {
    G4int bin=0;
    bin = G4int( (theNbin[l]-1)*(costh+1)/2. );
    if(bin == theNbin[l]-1) bin--;
    return bin;
  }
  
  inline G4double Interpolate(G4int bin, G4int l, G4double y1, G4double y2, G4double x)
  {
    G4double slope = 0, off = 0, x2=0, x1mx2;
    G4int half = (theNbin[l]-1)/2;
//    x1 = (bin-half)/G4double(half);
    x2 = (bin+1-half)/G4double(half);
    x1mx2 = 1./G4double( (theNbin[l]-1)/2 );
//    slope = (y2-y1)/(x2-x1);
    slope = (y2-y1)/x1mx2;
    off = y2-x2*slope;
    return x*slope+off;
  }
  
  const G4double ** value;
  const G4double ** integral;
  G4int theNbin[31];
  static const G4double l0[201];
  static const G4double i0[201];
  static const G4double l1[401];
  static const G4double i1[401];
  static const G4double l2[601];
  static const G4double i2[601];
  static const G4double l3[801];
  static const G4double i3[801];
  static const G4double l4[1001];
  static const G4double i4[1001];
  static const G4double l5[1201];
  static const G4double i5[1201];
  static const G4double l6[1401];
  static const G4double i6[1401];
  static const G4double l7[1601];
  static const G4double i7[1601];
  static const G4double l8[1801];
  static const G4double i8[1801];
  static const G4double l9[2001];
  static const G4double i9[2001];
  static const G4double l10[2201];
  static const G4double i10[2201];
  static const G4double l11[2401];
  static const G4double i11[2401];
  static const G4double l12[2601];
  static const G4double i12[2601];
  static const G4double l13[2801];
  static const G4double i13[2801];
  static const G4double l14[3001];
  static const G4double i14[3001];
  static const G4double l15[3201];
  static const G4double i15[3201];
  static const G4double l16[3401];
  static const G4double i16[3401];
  static const G4double l17[3601];
  static const G4double i17[3601];
  static const G4double l18[3801];
  static const G4double i18[3801];
  static const G4double l19[4001];
  static const G4double i19[4001];
  static const G4double l20[4201];
  static const G4double i20[4201];
  static const G4double l21[4401];
  static const G4double i21[4401];
  static const G4double l22[4601];
  static const G4double i22[4601];
  static const G4double l23[4801];
  static const G4double i23[4801];
  static const G4double l24[5001];
  static const G4double i24[5001];
  static const G4double l25[5201];
  static const G4double i25[5201];
  static const G4double l26[5401];
  static const G4double i26[5401];
  static const G4double l27[5601];
  static const G4double i27[5601];
  static const G4double l28[5801];
  static const G4double i28[5801];
  static const G4double l29[6001];
  static const G4double i29[6001];
  static const G4double l30[6201];
  static const G4double i30[6201];
};
#endif
