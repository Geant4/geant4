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
// By JPW, working, but to be cleaned up. @@@
// 22 Dec 2006 - DHW added isotope dependence
// G.Folger, 25-Nov-2009: extend to 100TeV, using a constant above 20GeV
//

#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadTmpUtil.hh"
#include "globals.hh"

G4NeutronInelasticCrossSection::G4NeutronInelasticCrossSection()
 : G4VCrossSectionDataSet("Wellisch-Laidlaw")
{}

G4NeutronInelasticCrossSection::~G4NeutronInelasticCrossSection()
{}

G4double G4NeutronInelasticCrossSection::
GetCrossSection(const G4DynamicParticle* aPart, 
                const G4Element* anEle, G4double /*aTemperature*/)
{
  G4int nIso = anEle->GetNumberOfIsotopes();
  G4double KE = aPart->GetKineticEnergy();
  G4double cross_section = 0;
  
  if (nIso) {
    G4double psig;
    G4IsotopeVector* isoVector = anEle->GetIsotopeVector();
    G4double* abundVector = anEle->GetRelativeAbundanceVector();
    G4int ZZ;
    G4int AA;

    for (G4int i = 0; i < nIso; i++) {
      ZZ = (*isoVector)[i]->GetZ();
      AA = (*isoVector)[i]->GetN();
      psig = GetCrossSection(KE, AA, ZZ);
      cross_section += psig*abundVector[i];
    }
  
  } else {
    G4int ZZ = G4lrint(anEle->GetZ());
    G4int AA = G4lrint(anEle->GetN());
    cross_section = GetCrossSection(KE, AA, ZZ);
  }

  return cross_section;
}

   
G4double G4NeutronInelasticCrossSection::
GetCrossSection(G4double anEnergy, G4int AA, G4int ZZ)
{
  G4double atomicNumber = G4double(AA);
  G4double nOfProtons = G4double(ZZ);

  if (anEnergy > 19.9*GeV ) 
  { // constant cross section above ~20GeV.
    return  GetCrossSection(19.8*GeV, AA, ZZ);
  } 
  G4double kineticEnergy = std::log10(DBL_MIN/MeV);
  if (anEnergy > DBL_MIN/MeV) kineticEnergy = std::log10(anEnergy/MeV);
  G4double nOfNeutrons = atomicNumber-nOfProtons;
  const G4double p1=1.3773;
  const G4double p2=1.+10./atomicNumber-0.0006*atomicNumber;
  const G4double p3=0.6+13./atomicNumber-0.0005*atomicNumber;
  const G4double p4=7.2449-0.018242*atomicNumber;
  const G4double p5=1.64-1.8/atomicNumber-0.0005*atomicNumber;
  const G4double p6=1.+200./atomicNumber+0.02*atomicNumber;
  const G4double p7=(atomicNumber-70.)*(atomicNumber-200.)/11000.;

  G4double logN = 1.0;
  if (nOfNeutrons > 1.5) logN = std::log(nOfNeutrons);      
  G4double part1 = pi*(p1*p1)*logN;
  G4double part2 = 1.+ std::pow(atomicNumber, 1./3.) 
                     - p2*(1.-1./std::pow(atomicNumber, 1./3.));

  G4double firstexp = -p4*(kineticEnergy-p5);
  G4double first=1.+std::exp(firstexp);
  G4double corr = 1.+p3*(1.-1./first); 

  G4double secondexp = -p6*(kineticEnergy-p7);
  G4double second=1.+std::exp(secondexp);
  G4double corr2 =1./second;

  G4double xsec = corr*corr2*part1*part2*10.*millibarn;
  return xsec;
}
