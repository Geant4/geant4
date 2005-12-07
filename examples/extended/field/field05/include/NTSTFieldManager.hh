#ifndef NTSTFieldMgr_hh
#define NTSTFieldMgr_hh

#include "G4FieldManager.hh"
class G4Track;

// Choose between two different sets of values for the FieldManager
//
  
class NTSTFieldManager : public G4FieldManager { 
     public:
        NTSTFieldManager(G4Field        *commonField,
			 G4FieldManager *pfmPrimary, 
			 G4FieldManager *pfmAlternate,
			 G4double       thresholdKineticEnergy);
        virtual ~NTSTFieldManager() {;}

        virtual void  ConfigureForTrack( const G4Track * ); 

     protected:	
        // inline (to be)
        const G4FieldManager& CopyValuesAndChordFinder( const G4FieldManager& newFieldMgr ); 
          // Copy the accuracy parameters and pointer to ChordFinder; 
          //  do not copy the field pointer -- for now!
     private:
        // Depending on Configuration, copy values from the following FM:
        G4FieldManager *pPrimaryFieldMgr;
        G4FieldManager *pAlternativeFieldMgr;

        G4double fThresholdKineticEnergy;
	// Threshold Kinetic energy:  if E < E_thr use alternative

        G4FieldManager *fpCurrentFieldMgr;   
            // Cached pointer to current one -- to avoid extra copies?
};

#endif  
