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
// $Id: G4RToEConvForGamma.cc 70745 2013-06-05 10:54:00Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4RToEConvForGamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

G4RToEConvForGamma::G4RToEConvForGamma()  
  : G4VRangeToEnergyConverter(),
    Z(-1),  
    s200keV(0.), s1keV(0.),
    tmin(0.),    tlow(0.), 
    smin(0.),    slow(0.),
    cmin(0.),    clow(0.), chigh(0.)
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  if (theParticle ==0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4RToEConvForGamma::G4RToEConvForGamma() ";
      G4cout << " Gamma is not defined !!" << G4endl;
    }
#endif
  } 
}

G4RToEConvForGamma::~G4RToEConvForGamma()
{ 
}


// ***********************************************************************
// ******************* BuildAbsorptionLengthVector ***********************
// ***********************************************************************
void G4RToEConvForGamma::BuildAbsorptionLengthVector(
                            const G4Material* aMaterial,
                            G4RangeVector* absorptionLengthVector )
{
  // fill the absorption length vector for this material
  // absorption length is defined here as
  //
  //    absorption length = 5./ macroscopic absorption cross section
  //
  const G4CrossSectionTable* aCrossSectionTable = (G4CrossSectionTable*)(theLossTable);
  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  //  fill absorption length vector
  G4int NumEl = aMaterial->GetNumberOfElements();
  G4double absorptionLengthMax = 0.0;
  for (size_t ibin=0; ibin<size_t(TotBin); ibin++) {
    G4double SIGMA = 0. ;
    for (size_t iel=0; iel<size_t(NumEl); iel++) {
      G4int IndEl = (*elementVector)[iel]->GetIndex();
      SIGMA +=  atomicNumDensityVector[iel]*
	           (*((*aCrossSectionTable)[IndEl]))[ibin];
    }
    //  absorption length=5./SIGMA
    absorptionLengthVector->PutValue(ibin, 5./SIGMA);
    if (absorptionLengthMax < 5./SIGMA ) absorptionLengthMax = 5./SIGMA;
  }
}



// ***********************************************************************
// ********************** ComputeCrossSection ****************************
// ***********************************************************************
G4double G4RToEConvForGamma::ComputeCrossSection(G4double AtomicNumber,
						 G4double KineticEnergy) 
{
  //  Compute the "absorption" cross section of the photon "absorption"
  //  cross section means here the sum of the cross sections of the
  //  pair production, Compton scattering and photoelectric processes
  const  G4double t1keV = 1.*keV;
  const  G4double t200keV = 200.*keV;
  const  G4double t100MeV = 100.*MeV;

  //  compute Z dependent quantities in the case of a new AtomicNumber
  if(std::abs(AtomicNumber-Z)>0.1)  {
    Z = AtomicNumber;
    G4double Zsquare = Z*Z;
    G4double Zlog = std::log(Z);
    G4double Zlogsquare = Zlog*Zlog;

    s200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;
    tmin = (0.552+218.5/Z+557.17/Zsquare)*MeV;
    smin = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*std::exp(1.5*Zlog);
    cmin=std::log(s200keV/smin)/(std::log(tmin/t200keV)*std::log(tmin/t200keV));
    tlow = 0.2*std::exp(-7.355/std::sqrt(Z))*MeV;
    slow = s200keV*std::exp(0.042*Z*std::log(t200keV/tlow)*std::log(t200keV/tlow));
    s1keV = 300.*Zsquare;
    clow =std::log(s1keV/slow)/std::log(tlow/t1keV);

    chigh=(7.55e-5-0.0542e-5*Z)*Zsquare*Z/std::log(t100MeV/tmin);
  }

  //  calculate the cross section (using an approximate empirical formula)
  G4double xs;
  if ( KineticEnergy<tlow ) {
    if(KineticEnergy<t1keV) xs = slow*std::exp(clow*std::log(tlow/t1keV));
    else                    xs = slow*std::exp(clow*std::log(tlow/KineticEnergy));

  } else if ( KineticEnergy<t200keV ) {
    xs = s200keV
         * std::exp(0.042*Z*std::log(t200keV/KineticEnergy)*std::log(t200keV/KineticEnergy));

  } else if( KineticEnergy<tmin ){
    xs = smin
         * std::exp(cmin*std::log(tmin/KineticEnergy)*std::log(tmin/KineticEnergy));

  } else {
    xs = smin + chigh*std::log(KineticEnergy/tmin);

  }
  return xs * barn;
}

