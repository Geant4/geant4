#ifndef MLNIELFunction_h
#define MLNIELFunction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4DataInterpolation.hh"

#include "MLNIELType.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLNIELFunction
{
public:
  MLNIELFunction();
  ~MLNIELFunction();

public:

  G4double GetNielFunction(G4int, NIELType, G4double);
//
//
// Density is for the moment hardwired for silicon in units of g/cm3.
//
  G4double GetDensity (NIELType nielType) {return 2.3290;};

private:  // NIEL functions data
  static  G4double RoseProtonEnergy[908] ;
  static  G4double RoseProton[908] ;
  static  G4double RoseNeutronEnergy[1381] ;
  static  G4double RoseNeutron[1381] ;
  static  G4double RosePionEnergy[900] ;
  static  G4double RosePion[900] ;
  static  G4double RoseElectronEnergy[15] ;
  static  G4double RoseElectron[15] ;

  static  G4double JPLProtonEnergy[28] ;
  static  G4double JPLProton[28] ;

  G4DataInterpolation *RoseProtonNiel;
  G4DataInterpolation *RoseNeutronNiel;
  G4DataInterpolation *RoseElectronNiel;
  G4DataInterpolation *RosePionNiel;

  G4DataInterpolation *JPLProtonNiel;

};
////////////////////////////////////////////////////////////////////////////////
#endif
