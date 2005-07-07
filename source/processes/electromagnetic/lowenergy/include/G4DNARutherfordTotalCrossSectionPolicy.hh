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
// $Id: G4DNARutherfordTotalCrossSectionPolicy.hh,v 1.2 2005-07-07 16:37:47 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNARUTHERFORDTOTALCROSSSECTIONPOLICY_HH
 #define  G4DNARUTHERFORDTOTALCROSSSECTIONPOLICY_HH 1
 
 #include "globals.hh"
 
 // EnergyLimitsPolicy must provide:
 //  - [public] static const double lowEnergyLimit
 //  - [public] static const double zeroBelowLowEnergyLimit
 //  - [public] static const double highEnergyLimit
 //  - [public] static const double zeroAboveLowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNARutherfordTotalCrossSectionPolicy
 {
  protected:
                                        G4DNARutherfordTotalCrossSectionPolicy() {}
                                       ~G4DNARutherfordTotalCrossSectionPolicy() {}
 
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
   G4double                             TotalCrossSection(G4double k, G4int z) const;
   void                                 BuildTotalCrossSection(void) const {}

  private:
   G4double                             RutherfordTotalCrossSection(G4double k, G4int z) const;
   G4double                             ScreeningFactor(G4double k, G4int z) const;

   // Hides default constructor and assignment operator as private 
                                        G4DNARutherfordTotalCrossSectionPolicy(const G4DNARutherfordTotalCrossSectionPolicy & copy);
   G4DNARutherfordTotalCrossSectionPolicy & operator=(const G4DNARutherfordTotalCrossSectionPolicy & right);
 };
 
 #include "G4DNARutherfordTotalCrossSectionPolicy.icc"
#endif /* G4DNARUTHERFORDTOTALCROSSSECTIONPOLICY_HH */

