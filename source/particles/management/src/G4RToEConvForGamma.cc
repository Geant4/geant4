//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RToEConvForGamma.cc,v 1.2 2002-12-16 11:15:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/iomanip"
#include "g4std/strstream"

G4RToEConvForGamma::G4RToEConvForGamma() : G4VRangeToEnergyConverter()
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
                            G4double       ,     
                            G4double       ,
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
    G4double lowEdgeEnergy = absorptionLengthVector->GetLowEdgeEnergy(ibin);
    
    G4double SIGMA = 0. ;
    
    for (size_t iel=0; iel<size_t(NumEl); iel++) {
      G4bool isOut;
      G4int IndEl = (*elementVector)[iel]->GetIndex();
      SIGMA +=  atomicNumDensityVector[iel]*
	           (*aCrossSectionTable)[IndEl]->GetValue(lowEdgeEnergy,isOut);
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
						 G4double KineticEnergy) const
{
  //  Compute the "absorption" cross section of the photon "absorption"
  //  cross section means here the sum of the cross sections of the
  //  pair production, Compton scattering and photoelectric processes
  static G4double Z;  
  const  G4double t1keV = 1.*keV;
  const  G4double t200keV = 200.*keV;
  const  G4double t100MeV = 100.*MeV;

  static G4double s200keV, s1keV;
  static G4double tmin, tlow; 
  static G4double smin, slow;
  static G4double cmin, clow, chigh;
  //  compute Z dependent quantities in the case of a new AtomicNumber
  if(abs(AtomicNumber-Z)>0.1)  {
    Z = AtomicNumber;
    G4double Zsquare = Z*Z;
    G4double Zlog = log(Z);
    G4double Zlogsquare = Zlog*Zlog;

    s200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;
    tmin = (0.552+218.5/Z+557.17/Zsquare)*MeV;
    smin = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*exp(1.5*Zlog);
    cmin=log(s200keV/smin)/(log(tmin/t200keV)*log(tmin/t200keV));
    tlow = 0.2*exp(-7.355/sqrt(Z))*MeV;
    slow = s200keV*exp(0.042*Z*log(t200keV/tlow)*log(t200keV/tlow));
    s1keV = 300.*Zsquare;
    clow =log(s1keV/slow)/log(tlow/t1keV);

    chigh=(7.55e-5-0.0542e-5*Z)*Zsquare*Z/log(t100MeV/tmin);
  }

  //  calculate the cross section (using an approximate empirical formula)
  G4double s;
  if ( KineticEnergy<tlow ) {
    if(KineticEnergy<t1keV) s = slow*exp(clow*log(tlow/t1keV));
    else                    s = slow*exp(clow*log(tlow/KineticEnergy));

  } else if ( KineticEnergy<t200keV ) {
    s = s200keV
         * exp(0.042*Z*log(t200keV/KineticEnergy)*log(t200keV/KineticEnergy));

  } else if( KineticEnergy<tmin ){
    s = smin
         * exp(cmin*log(tmin/KineticEnergy)*log(tmin/KineticEnergy));

  } else {
    s = smin + chigh*log(KineticEnergy/tmin);

  }
  return s * barn;
}

