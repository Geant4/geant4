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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITTRACKINGINTERACTIVITY_HH
#define G4ITTRACKINGINTERACTIVITY_HH

#include "G4String.hh"

class G4Track;
class G4Step;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4VITSteppingVerbose;

class G4ITTrackingInteractivity
{
protected:
  int fVerboseLevel;

public:
  G4ITTrackingInteractivity(G4VITSteppingVerbose* verbose = 0);

  virtual ~G4ITTrackingInteractivity()
  {
    ;
  }

  virtual void Initialize()
  {
    ;
  }

  virtual void StartTracking(G4Track* /*track*/)
  {

  }

  virtual void AppendStep(G4Track* /*track*/,
                          G4Step* /*step*/)
  {
    ;
  }

  virtual void EndTracking(G4Track* /*track*/)
  {

  }

  virtual void Finalize()
  {
    ;
  }

  virtual void TrackBanner(G4Track* /*track*/,
                           const G4String& message = "");

  inline void SetVerbose(int flag)
  {
    fVerboseLevel = flag;
  }

  inline G4int GetVerboseLevel() const
  {
    return fVerboseLevel;
  }

  void SetSteppingVerboseLevel(G4int level);
  G4int GetSteppingVerboseLevel() const;

  inline G4VITSteppingVerbose* GetSteppingVerbose()
  {
    return fpVerbose;
  }

  inline void SetSteppingVerbose(G4VITSteppingVerbose* verbose)
  {
    fpVerbose = verbose;
  }

private:
  G4VITSteppingVerbose* fpVerbose;
};

#endif // G4ITTRACKINGINTERACTIVITY_HH
