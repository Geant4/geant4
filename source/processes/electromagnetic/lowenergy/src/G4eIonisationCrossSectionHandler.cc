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
// $Id: G4eIonisationCrossSectionHandler.cc,v 1.8 2002-07-19 17:32:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eIonisationCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 25 Sept 2001
//
// Modifications: 
// 10 Oct 2001  M.G. Pia    Revision to improve code quality and consistency with design
// 19 Jul 2002   VI         Create composite data set for material
//
// -------------------------------------------------------------------

#include "G4eIonisationCrossSectionHandler.hh"
#include "G4VEnergySpectrum.hh"
#include "G4DataVector.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LinLogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"


G4eIonisationCrossSectionHandler::G4eIonisationCrossSectionHandler(
    const G4VEnergySpectrum* spec, G4VDataSetAlgorithm* alg, 
    G4double emin, G4double emax, G4int nbin)
 :  G4VCrossSectionHandler(),
    theParam(spec)
{
  G4VCrossSectionHandler::Initialise(alg, emin, emax, nbin);
  interp = new G4LinLogLogInterpolation();
}


G4eIonisationCrossSectionHandler::~G4eIonisationCrossSectionHandler() 
{
  delete interp;
}


G4std::vector<G4VEMDataSet*>* G4eIonisationCrossSectionHandler::BuildCrossSectionsForMaterials(
											       const G4DataVector& energyVector, 
											       const G4DataVector* energyCuts)
{
  G4int verbose = 0;
  G4std::vector<G4VEMDataSet*>* set = new G4std::vector<G4VEMDataSet*>;

  G4DataVector* energies;
  G4DataVector* cs;
  G4int nOfBins = energyVector.size();

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4VCrossSectionHandler::G4VCrossSectionHandler - no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material = (*materialTable)[m];
    const G4ElementVector* elementVector = material->GetElementVector();
    //const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
    const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();
    G4int nElements = material->GetNumberOfElements();

    if(verbose > 0) {
      G4cout << "eIonisation CS for " << m << "th material " 
             << material->GetName() 
             << "  eEl= " << nElements << G4endl; 
    }

    G4double tcut  = (*energyCuts)[m];
    
    G4VDataSetAlgorithm* algo = interp->Clone();
    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; i++) {
 
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4int nShells = NumberOfComponents(Z);
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nOfBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        G4double value = 0.0;

        if(e > tcut) {
          for (G4int n=0; n<nShells; n++) {
            G4double cross = FindValue(Z, e, n);
            G4double p = theParam->Probability(Z, tcut, e, e, n);  
            value += cross * p * density;
	    
            if(verbose>0 && m == 0 && e>=1. && e<=0.) {	    
              G4cout << "G4eIonCrossSH: e(MeV)= " << e/MeV
                     << " n= " << n
                     << " cross= " << cross
                     << " p= " << p
                     << " value= " << value
                     << " tcut(MeV)= " << tcut/MeV
                     << " rho= " << density
                     << " Z= " << Z
                     << G4endl;
	    }
	    
	  }
	}
        cs->push_back(value);
      } 
      G4VDataSetAlgorithm* algo = interp->Clone();
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,algo,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }
 
  return set;
}


