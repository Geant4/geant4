////////////////////////////////////////////////////////////////////////////////
//
#include "MLNIELFunction.hh"
#include "MLNIELData.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLNIELFunction::MLNIELFunction ()
{
  RoseProtonNiel   = new G4DataInterpolation(
                         RoseProtonEnergy, RoseProton, 908, 1e30, 1e30);
  RoseNeutronNiel  = new G4DataInterpolation(
                         RoseNeutronEnergy, RoseNeutron, 1381, 1e30, 1e30);
  RoseElectronNiel = new G4DataInterpolation(
                         RoseElectronEnergy, RoseElectron, 15, 1e30, 1e30);
  RosePionNiel     = new G4DataInterpolation(
                         RosePionEnergy, RosePion, 900, 1e30, 1e30);

  JPLProtonNiel    = new G4DataInterpolation(
                         JPLProtonEnergy, JPLProton, 28, 1e30, 1e30);
}
////////////////////////////////////////////////////////////////////////////////
//
MLNIELFunction::~MLNIELFunction ()
{ 
  delete RoseProtonNiel;
  delete RoseNeutronNiel;
  delete RosePionNiel;
  delete RoseElectronNiel;
  
  delete JPLProtonNiel;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double MLNIELFunction::GetNielFunction (G4int par, NIELType type,
  G4double ener)
{
//
//
// If the energy is zero, return zero NIEL.
//
  if (ener == 0.) return 0.;
//
//
// Initialise variables.  Units are MeV/(mg/cm2)  95*2.144e-5=2.0368e-6
//
  const G4double unitConv = 2.0368e-6;

  switch (par) {
  case 0 : // proton
    if (type == JPL) {
      if (ener <= JPLProtonEnergy[0] ) {
        return JPLProton[0]*1e-6;
      } else if ( ener >= JPLProtonEnergy[27] ) {
        return JPLProton[27]*1e-6;
      } else {
        return JPLProtonNiel->CubicSplineInterpolation(ener) * 1e-6 ; // original data in eV
      }
    } else {
      if (ener <= RoseProtonEnergy[0] ) {
        return RoseProton[0] * unitConv;
      } else if ( ener >= RoseProtonEnergy[907] ) {
        return RoseProton[907] * unitConv;
      } else {
        return RoseProtonNiel->CubicSplineInterpolation(ener) * unitConv;
      }
      //
    }
    break;
//
//
// Note that the remainder of this member functionshould not be reached if type
// == JPL, since the calling routine has an 'if' statement to calculate NIELs
// for neutrons, pions, and electrons if type == CERN.
//
  case 1 : // neutron
    if (ener <= RoseNeutronEnergy[0] ) {
      return RoseNeutron[0] * unitConv;
    } else if ( ener >= RoseNeutronEnergy[1380] ) {
      return RoseProton[1380] * unitConv;
    } else {
      return RoseNeutronNiel->CubicSplineInterpolation(ener) * unitConv;
    }
    break;
  case 2: // pion
    if (ener <= RosePionEnergy[0] ) {
      return RosePion[0] * unitConv;
    } else if ( ener >= RosePionEnergy[899] ) {
      return RosePion[899] * unitConv;
    } else {
      return RosePionNiel->CubicSplineInterpolation(ener) * unitConv;
    }
    break;
  case 3: // electron
    if (ener <= RoseElectronEnergy[0] ) {
      return RoseElectron[0] * unitConv;
    } else if ( ener >= RoseElectronEnergy[14] ) {
      return RoseElectron[14] * unitConv;
    } else {
      return RoseElectronNiel->CubicSplineInterpolation(ener) * unitConv;
    }
    break;
  }
//
//
// For satefy, the function returns -1.0e38 if end is reached.
//
  return -1.0e-38;
}
////////////////////////////////////////////////////////////////////////////////
