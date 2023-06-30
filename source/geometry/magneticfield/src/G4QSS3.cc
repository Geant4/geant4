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
// QSS3 implementation
//
// Author: Lucio Santi - 2020-2021.
// --------------------------------------------------------------------

#include "G4QSS3.hh"

G4QSS3::G4QSS3(QSS_simulator sim)
  : simulator(sim)
{
}

void G4QSS3::recompute_next_times(G4int *inf, G4double t)
{
  G4int i;
  G4double *x = simulator->x;
  G4double *q = simulator->q;
  G4double *lqu = simulator->lqu;
  G4double *time = simulator->nextStateTime;
  G4double coeff[4];

  for(i = 0; i < 3; i++)
  {
    const G4int var = inf[i],
                cf0 = 4*var,
                cf1 = cf0 + 1,
                cf2 = cf1 + 1,
                cf3 = cf2 + 1;

    if(std::fabs(q[cf0] - x[cf0]) >= lqu[var])
    {
      time[var] = t;
    }
    else
    {
      coeff[0] = q[cf0] + lqu[var] - x[cf0];
      coeff[1] = q[cf1] - x[cf1];
      coeff[2] = q[cf2] - x[cf2];
      coeff[3] = -x[cf3];
      time[var] = t + min_pos_root_3_alt(coeff, q[cf0] - lqu[var] - x[cf0]);
    }
  }
}
