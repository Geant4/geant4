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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLCascadeAction.hh
 * \brief Class containing default actions to be performed at intermediate cascade steps
 *
 * \date 22nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLCASCADEACTION_HH
#define G4INCLCASCADEACTION_HH 1

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  class CascadeAction {
    // class INCL must be a friend because it needs to call private methods
    friend class INCL;

    public:
    CascadeAction();
    virtual ~CascadeAction();

    virtual void beforeRunUserAction(Config const *) {}
    virtual void beforeCascadeUserAction(IPropagationModel *) {}
    virtual void beforePropagationUserAction(IPropagationModel *) {}
    virtual void beforeAvatarUserAction(IAvatar *, Nucleus *) {}
    virtual void afterAvatarUserAction(IAvatar *, Nucleus *, FinalState *) {}
    virtual void afterPropagationUserAction(IPropagationModel *, IAvatar *) {}
    virtual void afterCascadeUserAction(Nucleus *) {}
    virtual void afterRunUserAction() {}

    private:
    // These four methods should be private because the user must not be
    // allowed to override them
    void beforeRunAction(Config const *config);
    void beforeCascadeAction(IPropagationModel *);
    void beforePropagationAction(IPropagationModel *pm);
    void beforeAvatarAction(IAvatar *a, Nucleus *n);
    void afterAvatarAction(IAvatar *a, Nucleus *n, FinalState *fs);
    void afterPropagationAction(IPropagationModel *pm, IAvatar *avatar);
    void afterCascadeAction(Nucleus *);
    void afterRunAction();

    void beforeRunDefaultAction(Config const *config);
    void beforeCascadeDefaultAction(IPropagationModel *pm);
    void beforePropagationDefaultAction(IPropagationModel *pm);
    void beforeAvatarDefaultAction(IAvatar *a, Nucleus *n);
    void afterAvatarDefaultAction(IAvatar *a, Nucleus *n, FinalState *fs);
    void afterPropagationDefaultAction(IPropagationModel *pm, IAvatar *avatar);
    void afterCascadeDefaultAction(Nucleus *);
    void afterRunDefaultAction();

    private: // data members
    long stepCounter;
  };

}
#endif // G4INCLCASCADEACTION_HH
