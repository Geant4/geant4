// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronInelasticProcess.cc,v 1.1 1999-01-07 16:11:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Inelastic Process Class
 // J.L. Chuma, TRIUMF, 24-Mar-1997
 // Last modified: 27-Mar-1997
 // J.P. Wellisch: Bug hunting, 23-Apr-97
 // Modified by J.L.Chuma 8-Jul-97 to eliminate possible division by zero for sigma
//
// 14-APR-98 F.W.Jones: variant G4HadronInelastic process for
// G4CrossSectionDataSet/DataStore class design.
//
// 17-JUN-98 F.W.Jones: removed extraneous code causing core dump.
//
 
#include "G4HadronInelasticProcess.hh"
 
 G4double G4HadronInelasticProcess::GetMeanFreePath(
  const G4Track &aTrack,
  G4double previousStepSize,
  G4ForceCondition *condition )
  {
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    if( aParticle->GetDefinition() != theParticle )
      G4Exception( this->GetProcessName()+
                   " called for "+
                   aParticle->GetDefinition()->GetParticleName() );
    G4Material *aMaterial = aTrack.GetMaterial();
    G4int nElements = aMaterial->GetNumberOfElements();
    
    // returns the mean free path in GEANT4 internal units
    
    const RWTPtrVector<G4Element> *theElementVector =
      aMaterial->GetElementVector();
    
    const G4double *theAtomicNumDensityVector =
      aMaterial->GetAtomicNumDensityVector();
    
    G4Element *anElement = (*theElementVector)[0];
    G4int j = anElement->GetIndex();
    
    // This apparently should not be here (not useful and dumps core)
    // FWJ 17-JUN-1998
    //    G4bool isOutRange;
    //    G4double xSection = (*((*thePhysicsTable)(j))).GetValue(
    //     aParticle->GetTotalMomentum()/GeV, isOutRange );
    
    G4double sigma = 0.0;
    for( G4int i=0; i<nElements; ++i )
    {
      G4double xSection =
        GetMicroscopicCrossSection( aParticle, (*theElementVector)[i] );
      sigma += theAtomicNumDensityVector[i] * xSection;
    }
    if( sigma > 0.0 )
      return 1.0/sigma;
    else
      return DBL_MAX;
  }
 
 void 
  G4HadronInelasticProcess::BuildThePhysicsTable()
  {
   if (!theCrossSectionDataStore) {
     //      G4Exception("G4HadronInelasticProcess::BuildThePhysicsTable: "
     //                  "no CrossSectionDataStore");
      return;
   }

   theCrossSectionDataStore->BuildPhysicsTable(*theParticle);

   //    G4int numberOfElements = G4Element::GetNumberOfElements();
   //    thePhysicsTable = new G4PhysicsTable( numberOfElements );
   //    
   //    // make a PhysicsVector for each element
   //    
   //    static const G4ElementTable *theElementTable = G4Element::GetElementTable();
   //    for( G4int i=0; i<numberOfElements; ++i )
   //      (*thePhysicsTable)(i) =
   //        theCrossSectionData.MakePhysicsVector( *this, *theParticle,
   //                                               (*theElementTable)[i] );
  }
 
 G4double G4HadronInelasticProcess::GetMicroscopicCrossSection(
  const G4DynamicParticle *aParticle,
  const G4Element *anElement)
  {
    // returns the microscopic cross section in GEANT4 internal units
    
   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronInelasticProcess::GetMicroscopicCrossSection:"
                  "no CrossSectionDataStore");
      return DBL_MIN;
   }
   return theCrossSectionDataStore->GetCrossSection(aParticle, anElement);

   //    G4bool isOutRange;
   //    G4int j = anElement->GetIndex();
   //    
   //    G4double s = (*((*thePhysicsTable)(j))).GetValue(
   //     aParticle->GetTotalMomentum()/GeV, isOutRange );
   //    return s;
  }
 
 /* end of file */
