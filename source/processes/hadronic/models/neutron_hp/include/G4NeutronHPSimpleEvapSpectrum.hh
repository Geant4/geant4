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
// $Id: G4NeutronHPSimpleEvapSpectrum.hh,v 1.11 2006-06-29 20:49:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPSimpleEvapSpectrum_h
#define G4NeutronHPSimpleEvapSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPSimpleEvapSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPSimpleEvapSpectrum()
  {
    expm1 = std::exp(-1.);
  }
  ~G4NeutronHPSimpleEvapSpectrum()
  {
  }
  
  inline void Init(std::ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile,eV);
    theThetaDist.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  inline G4double Sample(G4double anEnergy) 
  {
    G4double theta = theThetaDist.GetY(anEnergy)*eV;
    G4double random, cut, max, result;
    max = 10.*theta;
    do
    {
      random = G4UniformRand();
      result = -theta*std::log(random); 
      cut = G4UniformRand();
    }
    while(cut>result/max); 
    return result;
  }
  
  private:
  
  inline G4double Evapo(G4double anEnergy, G4double theta)
  {
    G4double result = (anEnergy*eV)*std::exp(-anEnergy*eV/theta);
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theThetaDist;
  
};

#endif
