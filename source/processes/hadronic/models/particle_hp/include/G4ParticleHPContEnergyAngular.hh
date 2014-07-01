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
// $Id: G4ParticleHPContEnergyAngular.hh,v 1.1 2013/02/20 17:34:48 arce Exp $
// GEANT4 tag $Name: GAMOS-04-01-00 $
//
// 080721 Add ClearHistories() method by T. Koi
//
#ifndef G4ParticleHPContEnergyAngular_h
#define G4ParticleHPContEnergyAngular_h 1

#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ParticleHPContAngularPar.hh"
#include "G4InterpolationManager.hh"
class G4ParticleDefinition;

// we will need one of these per product.

class G4ParticleHPContEnergyAngular : public G4VParticleHPEnergyAngular
{
  public:
  
  G4ParticleHPContEnergyAngular(G4ParticleDefinition* proj)
    : theProjectile(proj)
  {
    theAngular = 0;
    currentMeanEnergy = -2;
  }
  
  ~G4ParticleHPContEnergyAngular()
  {
    if(theAngular!=0) delete [] theAngular;
  }
  
  public:
  
  void Init(std::istream & aDataFile)
  {
    aDataFile >> theTargetCode >> theAngularRep >> theInterpolation >> nEnergy;
    theAngular = new G4ParticleHPContAngularPar[nEnergy];
    theManager.Init(aDataFile);
    if( getenv("G4PHPTESTALLISOT") )  {
      if( theManager.GetScheme(0) != HISTO ) {// check only the first one, as all have the same interpolation scheme
	  G4cerr << " G4PHPTESTALLISOT interpolation scheme = " << theManager.GetScheme(0) << G4endl;
      }
    }



    //>GAMOS
    /*    if( theManager.GetNRanges() != 1 ) {
      G4cerr << " NRANGES = " <<  theManager.GetNRanges() << G4endl; // use G4UIcommand...
      throw G4HadronicException(__FILE__, __LINE__, ("nRanges is different than 1 ");   
      } */
    for(G4int i=0; i<nEnergy; i++)
    {
      theAngular[i].Init(aDataFile, theProjectile);
      theAngular[i].SetInterpolation(theInterpolation);
      if( i != 0 ) {
	theAngular[i].PrepareTableInterpolation(&(theAngular[i-1]));
      } else {
	theAngular[i].PrepareTableInterpolation((G4ParticleHPContAngularPar*)0);
      }
    }

    if( getenv("G4PHPTESTALLISOT") ) {
      //GAMOS CHECK IF INPUT ENERGIES HAVE DISCRETE OUTPUT ENERGIES WITH DIFFERENT XS
      std::map<double,double> discreteEXS;
      for( G4int ii = 0; ii < nEnergy; ii++ ) {
	G4int nDiscE = theAngular[ii].GetNDiscreteEnergies();
	if( nDiscE ) {
	  G4ParticleHPList * AngularVal = theAngular[ii].GetAngDataList();     
	  if( discreteEXS.size() == 0 ) {
	    for( G4int jj= 0; jj < nDiscE; jj++ ) {
	      discreteEXS[AngularVal[jj].GetLabel()] = AngularVal[jj].GetValue(0);
	    }
	  } else {
	    for( G4int jj= 0; jj < nDiscE; jj++ ) {
	      G4double discE = AngularVal[jj].GetLabel();
	      G4double XS = AngularVal[jj].GetValue(0);
	      std::map<double,double>::const_iterator itee = discreteEXS.find(discE);
	      if( itee == discreteEXS.end()  // OUPTUT DISCRETE ENERGY NOT FOUND IN PREVIOUS INPUT ENERGY
		  || XS != (*itee).second ) { // OUPTUT DISCRETE ENERGY XS <> IN PREVIOUS INPUT ENERGY
		for( itee = discreteEXS.begin(); itee != discreteEXS.end(); itee++ ) {
		  G4cerr << " ENERGYDISCRETE " << (*itee).first << " " << (*itee).second << G4endl;
		}
		G4cerr << " !!ENERGYDISCRETE NOT EQUAL IN PREVIOUS " << discE <<  " " << XS << G4endl;
	      } 
	    }
	  }
	}
      }
    }

  }
G4double MeanEnergyOfThisInteraction();
G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
  
  private:
  
  G4double theTargetCode;
  G4int theAngularRep;
  G4int nEnergy;
  
  G4int theInterpolation;

  G4InterpolationManager theManager; // knows the interpolation between stores
  G4ParticleHPContAngularPar * theAngular;
  
  G4double currentMeanEnergy;

  G4ParticleDefinition* theProjectile;

   public:
      void ClearHistories(); 
  
};
#endif
