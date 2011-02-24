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
// $Id: G4CrossSectionPairGG.cc,v 1.6 2011-01-09 02:37:48 dennis Exp $
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

G4CrossSectionPairGG::G4CrossSectionPairGG(G4VCrossSectionDataSet* low,
                                           G4double Etransit)
 : G4VCrossSectionDataSet("G4CrossSectionPairGG"),
   theLowX(low),  ETransition(Etransit)		    
{
  theHighX=new G4GlauberGribovCrossSection();
  verboseLevel=0;
}

G4CrossSectionPairGG::~G4CrossSectionPairGG()
{
//   The cross section registry wil delete these
//    delete theLowX;
//    delete theHighX;
}

G4bool G4CrossSectionPairGG::IsIsoApplicable(const G4DynamicParticle* aParticle,
                      G4int ZZ, G4int AA)
{
    G4bool isApplicable(false);
    G4double Ekin=aParticle->GetKineticEnergy();
    if (Ekin < ETransition ) 
    {
      isApplicable = theLowX->IsIsoApplicable(aParticle,ZZ,AA);
    } else {
      isApplicable = theHighX->IsIsoApplicable(aParticle,ZZ,AA);
    }
    
    return isApplicable;    
}

G4double G4CrossSectionPairGG::GetZandACrossSection(const G4DynamicParticle* aParticle,
                              G4int ZZ, G4int AA,
                              G4double aTemperature)
{
    G4double Xsec(0.);
    std::vector<ParticleXScale>::iterator iter;
    iter=scale_factors.begin();
    G4ParticleDefinition * pDef=aParticle->GetDefinition();
    while ( iter !=scale_factors.end() && (*iter).first != pDef ) {++iter;}
    
    
    G4double Ekin=aParticle->GetKineticEnergy();
    if (Ekin < ETransition ) 
    {
      Xsec=theLowX->GetZandACrossSection(aParticle,ZZ,AA,aTemperature);
    } else {
      Xsec=theHighX->GetInelasticGlauberGribov(aParticle,ZZ,AA)
           * (*iter).second[ZZ];
	if ( verboseLevel > 2 )
	{  G4cout << " scaling .." << ZZ << " " << AA << " " <<
		(*iter).second[ZZ]<< " " <<theHighX->GetInelasticGlauberGribov(aParticle,ZZ,AA) << "  " 
		<< Xsec << G4endl;
	}	   
    }
    
    return Xsec;
}



void G4CrossSectionPairGG::BuildPhysicsTable(const G4ParticleDefinition& pDef)
{
    theLowX->BuildPhysicsTable(pDef);
    theHighX->BuildPhysicsTable(pDef);
    
    G4NistManager* NistMan = G4NistManager::Instance();
    G4ParticleDefinition * myDef=const_cast<G4ParticleDefinition*>(&pDef);
    std::vector<ParticleXScale>::iterator iter;
    iter=scale_factors.begin();
    while ( iter !=scale_factors.end() && (*iter).first != myDef ) {++iter;}

    //  new particle, initialise
    
    if ( iter == scale_factors.end() )
    {
       XS_factors factors (93);
       G4ThreeVector mom(0.0,0.0,1.0);
       G4DynamicParticle DynPart(myDef, mom, ETransition);  // last is kinetic Energy

       if (verboseLevel > 0) { 
	  G4cout << "G4CrossSectionPairGG::BuildPhysicsTable for particle "
	         << pDef.GetParticleName() << G4endl;
       }	   
       for (G4int aZ=1; aZ<93; ++aZ)
       {
          factors[aZ]=1.;   // default, to give reasonable value if only high is applicable
          G4int AA=G4lrint(NistMan->GetAtomicMassAmu(aZ));
          G4bool isApplicable = theLowX->IsIsoApplicable(&DynPart, aZ, AA) &&
                      theHighX->IsIsoApplicable(&DynPart, aZ, AA-aZ);
		   
	  if (isApplicable)
	  {
	     factors[aZ]=theLowX->GetZandACrossSection(&DynPart,aZ,AA,0) /
	              theHighX->GetInelasticGlauberGribov(&DynPart,aZ,AA);
                      
	  }  
	  if (verboseLevel > 0) { 
	     G4cout << "Z=" << aZ<< ",  A="<< AA << ", scale=" << factors[aZ];
	     if ( verboseLevel == 1) { G4cout  << G4endl; }
	     else {
		if (isApplicable) {
		   G4cout << ",  low / high " <<  theLowX->GetZandACrossSection(&DynPart,aZ,AA,0)
	        	 << "  " << theHighX->GetInelasticGlauberGribov(&DynPart,aZ,AA) << G4endl;
		} else { G4cout << ",   N/A" << G4endl; }	 
             }
	  }
       }
       ParticleXScale forPart(myDef,factors);
       scale_factors.push_back(forPart);
    }
}


void G4CrossSectionPairGG::DumpPhysicsTable(const G4ParticleDefinition&)
{
  G4cout << std::setw(24) << " " << " G4CrossSectionPairGG: "
         << theLowX->GetName() << " cross sections " << G4endl;
  G4cout << std::setw(27) << " " << "below " << ETransition/GeV
         << " GeV, Glauber-Gribov above " << G4endl;
}
