// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FTFModel.hh,v 1.4 2000/12/14 09:25:54 hpw Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// Class Description
// Final state production code for hadron inelastic scattering above 20 GeV
// based on the modeling ansatz used in FRITIOF.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with an object
// of G4TheoFSGenerator. 
// Class Description - End

#ifndef G4FTFModel_h
#define G4FTFModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------


#include "G4VPartonStringModel.hh"

class G4VSplitableHadron;
class G4ExcitedString;
#include "G4FTFParticipants.hh"

#include "G4ExcitedStringVector.hh"
#include "G4DiffractiveExcitation.hh"


class G4FTFModel : public G4VPartonStringModel
{

  public:
      G4FTFModel(G4double sigmaPt=800*MeV, G4double minExtraMass=300*MeV,G4double x0Mass=150*MeV);
      G4FTFModel(G4DiffractiveExcitation * anExcitation);
      G4FTFModel(const G4FTFModel &right);
      ~G4FTFModel();
      const G4FTFModel & operator=(const G4FTFModel &right);

      int operator==(const G4FTFModel &right) const;
      int operator!=(const G4FTFModel &right) const;

      void Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile);
      G4ExcitedStringVector * GetStrings();
      G4V3DNucleus * GetWoundedNucleus() const;


  protected:
  
  private:
       G4bool ExciteParticipants();
       G4ExcitedStringVector * BuildStrings();
  
  private:     
       
       G4FTFParticipants theParticipants;
       G4ReactionProduct theProjectile;
       
       G4DiffractiveExcitation * theExcitation;



};

inline 
G4V3DNucleus * G4FTFModel::GetWoundedNucleus() const
{
	return theParticipants.GetWoundedNucleus();
}

#endif


