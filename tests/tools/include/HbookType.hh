//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef HbookType_hh
#define HbookType_hh

#include "HbookManager.hh"

class HbookType
{
public:

  HbookType()
  {
    theHbookManager.Register(this);
  }

  virtual ~HbookType()
  {
    theHbookManager.UnRegister(this);
  }    

  // To define this const isn´t fair but avoids 
  // anoying compiler warnings
  int GetHbookID() const
  {
    return id;
  }

protected:

  int id;

private:

  void SetID(int v)
  {
    id=v;
  }

  friend class HbookManager;
};

#endif
