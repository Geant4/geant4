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

#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4HadTmpUtil.hh"
#include "G4ElNucleusSFcs.hh"
#include "G4ElectroNuclearCrossSection.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ElNucleusSFcs);


G4ElNucleusSFcs::G4ElNucleusSFcs():G4VCrossSectionDataSet(Default_Name())
{
  fCHIPScs = new G4ElectroNuclearCrossSection();  
}

G4ElNucleusSFcs::~G4ElNucleusSFcs()
{
  if( fCHIPScs != nullptr ) delete fCHIPScs;
}

void
G4ElNucleusSFcs::CrossSectionDescription( std::ostream& outFile) const
{
    outFile << "G4ElNucleusSFcs provides the inelastic\n"
    << "cross section for e- and e+ interactions with nuclei according to the structure function approach."
    << "all energies.\n";
}

G4bool G4ElNucleusSFcs::IsElementApplicable( const G4DynamicParticle* aPart, G4int Z, const G4Material*)
{

  G4double eKin   = aPart->GetKineticEnergy();
  G4double thTkin = ThresholdEnergy();
  
  return ( Z > 0 && Z < 120 && eKin > thTkin );
}


G4double G4ElNucleusSFcs::GetIsoCrossSection( const G4DynamicParticle* aPart, G4int ZZ, G4int AA,  const G4Isotope* ,
                              const G4Element* , const G4Material* )
{
  G4double xsc(0.), ratio(1.);

  xsc = fCHIPScs->GetElementCrossSection( aPart, ZZ, nullptr ); //, const G4Material*);

  ratio = GetRatio(ZZ,AA);

  if( ratio > 0. ) xsc /= ratio;
  
  return xsc;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate ratio CHIPS-CS/SF-CS, averaged between 200 MeV and 100 GeV

G4double G4ElNucleusSFcs::GetRatio( G4int Z, G4int A )
{
  G4double ratio(1.);

  if      ( Z == 1   && A == 1 )   return 1.51;
  else if ( Z == 1   && A == 2 )   return 0.33;
  else if ( Z == 1   && A == 3 )   return 0.27;
  else if ( Z == 2   && A == 4 )   return 1.81;
  else if ( Z == 6   && A == 12 )  return 2.26;
  else if ( Z == 7   && A == 14 )  return 2.47;
  else if ( Z == 8   && A == 16 )  return 2.61;
  else if ( Z == 13  && A == 27 )  return 2.57;
  else if ( Z == 14  && A == 28 )  return 2.49;
  else if ( Z == 18  && A == 40 )  return 2.72;
  else if ( Z == 22  && A == 48 )  return 2.71;
  else if ( Z == 26  && A == 56 )  return 2.79;
  else if ( Z == 29  && A == 64 )  return 2.78;
  else if ( Z == 32  && A == 73 )  return 2.87;
  else if ( Z == 42  && A == 96 )  return 3.02;
  else if ( Z == 46  && A == 106 ) return 3.02;
  else if ( Z == 47  && A == 108 ) return 2.99;
  else if ( Z == 48  && A == 112 ) return 3.00;
  else if ( Z == 74  && A == 184 ) return 3.44;
  else if ( Z == 79  && A == 200 ) return 3.49;
  else if ( Z == 82  && A == 207 ) return 3.48;
  else if ( Z == 92  && A == 238 ) return 3.88;
  else
  {
    G4int it(0), iMax(19);
    G4double zz = G4double(Z);

    for ( it = 0; it < iMax; ++it ) if ( zz <= fZZ[it] ) break;

    if     ( it == 0 )    return fRR[0];
    else if( it == iMax ) return fRR[iMax-1];
    else
    {
      G4double x1 = fZZ[it-1];
      G4double x2 = fZZ[it];
      G4double y1 = fRR[it-1];
      G4double y2 = fRR[it];
      
      if( x1 >= x2 ) return fRR[it];
      else
      {
        G4double angle = (y2-y1)/(x2-x1);
	ratio = y1 + ( zz - x1 )*angle;
      }
    } 
  }
  return ratio;
}

/////////////////////////////////////////

G4double G4ElNucleusSFcs::ThresholdEnergy()
{
  G4double thTkin = 134.9766*CLHEP::MeV; // Pi0 threshold for the nucleon
  
  thTkin *= 1.3; // nucleon motion
  
  return thTkin;
}

/////////////////////////////////////////

const G4double G4ElNucleusSFcs::fZZ[19] =
  {
    2., 6., 7., 8., 13., 14., 18., 22., 26., 29.,
    32., 42., 46., 47., 48., 74., 79., 82., 92.
  };

/////////////////////////////////////////

const G4double G4ElNucleusSFcs::fRR[19] =
  {
    1.81, 2.26, 2.47, 2.61, 2.57, 2.49, 2.72, 2.71, 2.79, 2.78,
    2.87, 3.02, 3.02, 2.99, 3., 3.44, 3.49, 3.48, 2.88
  };
