// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05MaxTimeCuts.hh,v 1.2 1999-06-09 10:33:06 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// ------------------------------------------------------------
//                  14 Aug. 1998  H.Kurashige
// ------------------------------------------------------------

#ifndef ExN05MaxTimeCuts_h
#define ExN05MaxTimeCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "ExN05SpecialCuts.hh"


class ExN05MaxTimeCuts : public ExN05SpecialCuts
{
  public:     

     ExN05MaxTimeCuts(const G4String& processName ="ExN05MaxTimeCuts" );

     virtual ~ExN05MaxTimeCuts();

     // PostStep GPIL
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );
            
			    
  private:
  
  // hide assignment operator as private 
      ExN05MaxTimeCuts(ExN05MaxTimeCuts&);
      ExN05MaxTimeCuts& operator=(const ExN05MaxTimeCuts& right);

};

#endif

