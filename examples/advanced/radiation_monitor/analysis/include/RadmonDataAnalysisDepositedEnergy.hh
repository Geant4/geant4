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
// File name:     RadmonDataAnalysisDepositedEnergy.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisDepositedEnergy.hh,v 1.3 2006/06/29 16:06:53 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   Energy deposit analysis
//

#ifndef   RADMONDATAANALYSISDEPOSITEDENERGY_HH
 #define  RADMONDATAANALYSISDEPOSITEDENERGY_HH
 
 // Include files
 #include "RadmonVDataAnalysisWithLabel.hh"
 
 class RadmonDataAnalysisDepositedEnergy : public RadmonVDataAnalysisWithLabel
 {
  public:
                                                RadmonDataAnalysisDepositedEnergy();
   inline virtual                              ~RadmonDataAnalysisDepositedEnergy();

   virtual RadmonVDataAnalysisWithLabel *       New(void) const;

   virtual G4String                             ObtainColumnsDeclaration(const G4String & prefix);
   virtual void                                 InitializeFromTuple(const G4String & prefix, const AIDA::ITuple * tuple);
   virtual void                                 StoreIntoTuple(RadmonHitsCollection * hitsCollection, AIDA::ITuple * tuple);

   virtual void                                 StoreIntoHit(G4Step * step, RadmonHit * hit);
  private:
  // Hidden constructors and operators
                                                RadmonDataAnalysisDepositedEnergy(const RadmonDataAnalysisDepositedEnergy & copy);
   RadmonDataAnalysisDepositedEnergy &          operator=(const RadmonDataAnalysisDepositedEnergy & copy);

  // Private attributes
   G4int                                        indexHitEnergyDeposit;
   G4int                                        indexTupleEnergyDeposit;
 };
 
 #include "RadmonDataAnalysisDepositedEnergy.icc"
#endif /* RADMONDATAANALYSISDEPOSITEDENERGY_HH */
