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
// $Id: G4CrossSectionPairGG.cc,v 1.1 2010-11-16 13:38:15 gunter Exp $
// $ GEANT4 tag $Name: not supported by cvs2svn $
//
//   Class G4CrossSectionPairGG
//
//     smoothly join two cross section sets by scaling the second at a given 
//       transition energy to match the first.
//
//  Author:  Gunter Folger
//           November 2009
//

#include "G4CrossSectionPairGG.hh"

#include "globals.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"

G4CrossSectionPairGG::G4CrossSectionPairGG(G4VCrossSectionDataSet * low,
//  	                    G4VCrossSectionDataSet * high,
			    G4double Etransit):
	theLowX(low),
//	theHighX(high),
	ETransition(Etransit)		    
{
    theHighX=new G4GlauberGribovCrossSection();
}

G4CrossSectionPairGG::~G4CrossSectionPairGG()
{
    delete theLowX;
    delete theHighX;
}


G4bool G4CrossSectionPairGG::IsZAApplicable(const G4DynamicParticle* particle,
                      G4double ZZ, G4double AA)
{
    G4bool isApplicable(false);
    G4double Ekin=particle->GetKineticEnergy();
    if (Ekin < ETransition ) 
    {
      isApplicable = theLowX->IsZAApplicable(particle,ZZ,AA);
    } else {
      isApplicable = theHighX->IsZAApplicable(particle,ZZ,AA);
    }
    
    return isApplicable;    
}

G4double G4CrossSectionPairGG::GetIsoZACrossSection(const G4DynamicParticle* aParticle,
                              G4double ZZ, G4double AA,
                              G4double aTemperature)
{
    G4double Xsec(0.);
    std::vector<ParticleXScale>::iterator iter;
    iter=scale_factors.begin();
    G4ParticleDefinition * pDef=aParticle->GetDefinition();
    while ( iter !=scale_factors.end() && (*iter).first != pDef ) {}
    
    
    G4double Ekin=aParticle->GetKineticEnergy();
    if (Ekin < ETransition ) 
    {
      Xsec=theLowX->GetIsoZACrossSection(aParticle,ZZ,AA,aTemperature);
    } else {
      Xsec=theHighX->GetInelasticGlauberGribov(aParticle,ZZ,AA) *
           (*iter).second[G4lrint(ZZ)];
	G4cout << " scaling .." << ZZ << " " << AA << " " <<
	(*iter).second[G4lrint(ZZ)]<< " " 
	<<theHighX->GetInelasticGlauberGribov(aParticle,ZZ,AA) << "  " 
	<< Xsec << G4endl;   
    }
    
    return Xsec;
}



void G4CrossSectionPairGG::BuildPhysicsTable(const G4ParticleDefinition& pDef)
{
    G4NistManager* NistMan = G4NistManager::Instance();
    G4ParticleDefinition * myDef=const_cast<G4ParticleDefinition*>(&pDef);
    std::vector<ParticleXScale>::iterator iter;
    iter=scale_factors.begin();
    while ( iter !=scale_factors.end() && *(*iter).first != pDef ) {}

    //  new particle, initialise
    
    if ( iter == scale_factors.end() )
    {
       XS_factors factors (93);
       G4ThreeVector mom(0.0,0.0,1.0);
       G4DynamicParticle DynPart(myDef, mom, ETransition);  // last is kinetic Energy

       for (G4int aZ=1; aZ<93; ++aZ)
       {
          factors[aZ]=1.;   // default, to give reasonable value if only high is applicable
          G4double AA=NistMan->GetAtomicMassAmu(aZ);
          G4bool isApplicable = theLowX->IsZAApplicable(&DynPart, G4double(aZ), AA) &&
                      theHighX->IsZAApplicable(&DynPart, G4double(aZ), AA);
		   
	  if (isApplicable)
	  {
	     factors[aZ]=theLowX->GetIsoZACrossSection(&DynPart,G4double(aZ),AA,0) /
	              theHighX->GetInelasticGlauberGribov(&DynPart,G4double(aZ),AA);
//	              theHighX->GetIsoZACrossSection(&DynPart,G4double(aZ),AA,0);
                      
	  }  
	G4cout << " xs scale : " << aZ << " " << AA << " " << factors[aZ]<< G4endl;
	G4cout << "			low / high " << theLowX->GetIsoZACrossSection(&DynPart,G4double(aZ),AA,0) 
	               << "  " << theHighX->GetInelasticGlauberGribov(&DynPart,G4double(aZ),AA) << G4endl;
       }
       ParticleXScale forPart(myDef,factors);
       scale_factors.push_back(forPart);
    }
}
void G4CrossSectionPairGG::DumpPhysicsTable(const G4ParticleDefinition&)
{
}
