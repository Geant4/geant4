//
// File name:     RadmonDataAnalysisDepositedEnergy.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisDepositedEnergy.hh,v 1.1 2005-11-25 01:56:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
