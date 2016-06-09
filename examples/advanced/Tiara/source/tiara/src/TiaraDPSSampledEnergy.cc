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
// $Id: TiaraDPSSampledEnergy.cc,v 1.9 2006/06/29 15:44:48 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#include "TiaraDPSSampledEnergy.hh"

#ifdef G4ANALYSIS_USE

#include "Randomize.hh"
#include "TiaraFileAcess.hh"


#include "AIDA/AIDA.h"


TiaraDPSSampledEnergy::TiaraDPSSampledEnergy(const G4String &eng,
					     G4double minEnergyCut,
					     const G4String &sourceTree,
					     const G4String &nameExt)
  :
 //  fTree( (checkFileIsReadable(sourceTree, "TiaraDPSSampledEnergy::TiaraDPSSampledEnergy"),
// 	 AIDA_createAnalysisFactory()
// 	->createTreeFactory()
// 	->create(sourceTree, "xml",true,false))  ),
  //fSampleDPS(dynamic_cast<AIDA::IDataPointSet *>(fTree->find(G4String("dps" + eng + nameExt)))),
  fMinEnergyCut(minEnergyCut),
  fMaxProb(0)
  // fMinE(fSampleDPS->point(0)->coordinate(0)->value()),
  //fMaxE(fSampleDPS->point(fSampleDPS->size()-1)->coordinate(0)->value())
{
  /////////////////////////////////////////
 
  AIDA::IAnalysisFactory* aFact = AIDA_createAnalysisFactory();

  AIDA::ITreeFactory *treeFact = aFact->createTreeFactory(); 
 
  fTree = treeFact -> create(sourceTree, "xml",true,false); 

  AIDA::IManagedObject* object = fTree->find("dps" + eng + nameExt);
  if(object) {
    AIDA::IDataPointSet* sample = object->cast("AIDA::IDataPointSet");
    fSampleDPS = sample;
  } 
  fMinE = fSampleDPS->point(0)->coordinate(0)->value();
  fMaxE = fSampleDPS->point(fSampleDPS->size()-1)->coordinate(0)->value();

  /////////////////////////////////////////

  for (G4int i=0; i < fSampleDPS->size(); i++) {
    AIDA::IDataPoint *p = fSampleDPS->point(i);
    G4double prob(p->coordinate(1)->value());
    fEnergy_Flux[p->coordinate(0)->value() * 10] = prob;
    if (prob > fMaxProb) {
      fMaxProb = prob;
    }
  }
  if (!(fMaxProb > 0)) {
    G4Exception("TiaraDPSSampledEnergy::TiaraDPSSampledEnergy: fMaxProb <= 0!");
  }
}

TiaraDPSSampledEnergy::TiaraDPSSampledEnergy(const TiaraDPSSampledEnergy& rhs)
  : TiaraVSourceEnergyGenerator()
{
  *this = rhs;
}

TiaraDPSSampledEnergy& TiaraDPSSampledEnergy::
operator=(const TiaraDPSSampledEnergy& rhs) {
  if (this!=&rhs) {
    TiaraDPSSampledEnergy &nonConstRhs = static_cast<TiaraDPSSampledEnergy &>(rhs);
    fTree = nonConstRhs.fTree;
    fSampleDPS = nonConstRhs.fSampleDPS;
    fEnergy_Flux = nonConstRhs.fEnergy_Flux;
    fMinEnergyCut = nonConstRhs.fMinEnergyCut;
    fMaxProb = nonConstRhs.fMaxProb;
    fMinE = nonConstRhs.fMinE;
    fMaxE = nonConstRhs.fMaxE;
  }
  return *this;
}

TiaraDPSSampledEnergy::~TiaraDPSSampledEnergy(){
}


TiaraVSourceEnergyGenerator *TiaraDPSSampledEnergy::Clone() const {
  return new TiaraDPSSampledEnergy(*this);
}

G4double TiaraDPSSampledEnergy::GetEnergy() {

  G4double e(0.0);
  while (!(e>fMinEnergyCut)) { 
    // a factor 10 is introduced to have an integer map 
    // for .5 values
    G4double r(10 * (fMinE + G4UniformRand() * (fMaxE - fMinE) ));
    G4int cLow(0);
    G4int cHigh(0);
    getBounds(cLow,cHigh,r);
    G4double dx10(r-cLow); // a facto 10 enters 

    G4double fLow(fEnergy_Flux[cLow]);
    G4double fHigh(fEnergy_Flux[cHigh]);

    //    if ( (!(fLow>0.0)) || (!(fHigh>0.0)) ) {
    //      G4cout << "TiaraDPSSampledEnergy::GetEnergy: WARNING: (!(fLow>0.0) || (!(fHigh>0.0) " << G4endl;
    //      G4cout << "r: " << r << G4endl;
    //      G4cout << "d: " << d << G4endl;
    //     G4cout << "cLow: " << cLow << G4endl;
    //      G4cout << "cHigh: " << cHigh << G4endl;
    //    }
    //    G4cout << "fLow: " << fLow << ", fHigh: " << fHigh << G4endl;

    G4double dy(fHigh - fLow);
    G4double m01(dy/(cHigh-cLow)); // a factor 1/10 enters
    G4double p((fLow + m01*dx10)/fMaxProb);
    //    G4cout << "p: " << p << G4endl;
    G4double discard(G4UniformRand());
    //    G4cout << "discard: " << discard << G4endl;
    if (discard < p) {
      e = r/10.;
      break;
    }
  }
  return e;
}

void TiaraDPSSampledEnergy::getBounds(G4int &cL, G4int &cH, G4double v) {
  std::map<int, double>::iterator itH = fEnergy_Flux.upper_bound(v);
  cH = itH->first;
  std::map<int, double>::iterator itL = --itH;
  cL = itL->first;
}

#endif







