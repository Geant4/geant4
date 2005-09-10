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
// $Id: G4DNATotalCrossSectionFromFilePolicy.hh,v 1.3 2005-09-10 09:34:33 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
 
   G4double                             TotalCrossSection(G4double k, G4int z) const;
   G4int                                RandomizePartialCrossSection(G4double k, G4int z);
   G4int                                NumberOfPartialCrossSections(void);
   void                                 BuildTotalCrossSection(void);

  private:
   void                                 Free(void);
  
   G4DNACrossSectionDataSet *           dataset;
   G4double *                           valuesBuffer;
   DataFilePolicy                       dataFilePolicy;

   // Hides default constructor and assignment operator as private 
                                        G4DNATotalCrossSectionFromFilePolicy(const G4DNATotalCrossSectionFromFilePolicy & copy);
   G4DNATotalCrossSectionFromFilePolicy & operator=(const G4DNATotalCrossSectionFromFilePolicy & right);
 };
 
 #include "G4DNATotalCrossSectionFromFilePolicy.icc"
#endif /* G4DNATOTALCROSSSECTIONFROMFILEPOLICY_HH */

