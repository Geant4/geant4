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
// $Id: G4ITBox.hh 100802 2016-11-02 14:55:27Z gcosmo $
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

#ifndef G4ITBox_h
#define G4ITBox_h

#include "G4IT.hh"

/**
 * G4ITBox behaves just like a stack for G4IT.
 * You can search for specific tracks.
 * Each G4IT knows to which G4ITBox it belongs and its corresponding node.
 * This makes the deletion of an element very fast.
 * The drawback is that a G4IT can only belong to one G4ITBox.
 * If you are not looking for this feature, please use std::list.
 */

class G4ITBox
{

public:
  G4ITBox();
  ~G4ITBox();

  void ResetStack();
  void Push(G4IT*);
  void Extract(G4IT*);

  /** The FindIT methods are used for check only.
   * Those methods are not effective due to the
   * linear search. It is better to use GetIT(track)
   * in order to retrieve the IT and GetIT(track)->GetBox()
   * in order to check which is the box pointer.
   */
  G4IT* FindIT(const G4Track&);
  const G4IT* FindIT(const G4Track&) const;
  void TransferTo(G4ITBox*);

  inline G4bool Empty() const;
  inline G4int GetNTrack() const;

  inline G4IT* GetFirstIT();
  inline G4IT* GetLastIT();
  inline const G4IT* GetFirstIT() const;
  inline const G4IT* GetLastIT() const;

  inline void SetNextBox(G4ITBox* box);
  inline G4ITBox* GetNextBox();
  inline const G4ITBox* GetNextBox() const;
  inline void SetPreviousBox(G4ITBox* box);
  inline G4ITBox* GetPreviousBox();
  inline const G4ITBox* GetPreviousBox() const;

private:
  const G4ITBox & operator=(const G4ITBox &right);
  G4int fNbIT;
  G4IT * fpFirstIT;
  G4IT * fpLastIT;

  G4ITBox* fpPreviousBox;
  G4ITBox* fpNextBox;
};

inline G4bool G4ITBox::Empty() const
{
  return (fNbIT <= 0);
}

inline G4int G4ITBox::GetNTrack() const
{
  return fNbIT;
}
inline G4IT* G4ITBox::GetFirstIT()
{
  return fpFirstIT;
}
inline G4IT* G4ITBox::GetLastIT()
{
  return fpLastIT;
}

inline const G4IT* G4ITBox::GetFirstIT() const
{
  return fpFirstIT;
}
inline const G4IT* G4ITBox::GetLastIT() const
{
  return fpLastIT;
}

inline void G4ITBox::SetNextBox(G4ITBox* box)
{
  fpNextBox = box;
}

inline G4ITBox* G4ITBox::GetNextBox()
{
  return fpNextBox;
}

inline void G4ITBox::SetPreviousBox(G4ITBox* box)
{
  fpPreviousBox = box;
}

inline G4ITBox* G4ITBox::GetPreviousBox()
{
  return fpPreviousBox;
}

#endif
