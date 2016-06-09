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
// $Id: G4DNATotalCrossSectionFromFilePolicy.hh,v 1.5 2006/06/29 19:35:27 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

#ifndef   G4DNATOTALCROSSSECTIONFROMFILEPOLICY_HH
#define  G4DNATOTALCROSSSECTIONFROMFILEPOLICY_HH 1
 
#include "G4DNACrossSectionDataSet.hh"
 
// IncomingParticlePolicy must provide:
//  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);
 
// DataFilePolicy must provide:
//  - [public] static const double lowEnergyLimit
//  - [public] static const double zeroBelowLowEnergyLimit
//  - [public] static const double highEnergyLimit
//  - [public] static const double zeroAboveLowEnergyLimit
//  - [public] static const double dataFileEnergyUnit
//  - [public] static const double dataFileCrossSectionUnit
//  - [public] static char const * const dataFileName
 
// InterpolationAlgorithmPolicy must inherit from [public] G4VDataSetAlgorithm

template <typename IncomingParticlePolicy, typename DataFilePolicy, typename InterpolationAlgorithmPolicy>
class G4DNATotalCrossSectionFromFilePolicy : public IncomingParticlePolicy
{
protected:
  G4DNATotalCrossSectionFromFilePolicy();
  ~G4DNATotalCrossSectionFromFilePolicy();
 
  G4double TotalCrossSection(G4double k, G4int z) const;
  G4int RandomizePartialCrossSection(G4double k, G4int z);
  G4int NumberOfPartialCrossSections(void);
  void  BuildTotalCrossSection(void);

private:
  void Free(void);
  
  G4DNACrossSectionDataSet* dataset;
  G4double* valuesBuffer;
  DataFilePolicy dataFilePolicy;

  // Hides default constructor and assignment operator as private 
  G4DNATotalCrossSectionFromFilePolicy(const G4DNATotalCrossSectionFromFilePolicy & copy);
  G4DNATotalCrossSectionFromFilePolicy & operator=(const G4DNATotalCrossSectionFromFilePolicy & right);
};
 
#include "G4DNATotalCrossSectionFromFilePolicy.icc"
#endif /* G4DNATOTALCROSSSECTIONFROMFILEPOLICY_HH */

