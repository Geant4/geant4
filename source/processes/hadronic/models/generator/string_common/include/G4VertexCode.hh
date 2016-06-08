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
//
// $Id: G4VertexCode.hh,v 1.5 2002/12/12 19:17:55 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
#ifndef G4VertexCode_h
#define G4VertexCode_h 1


class G4VertexCode
{

  public:
    G4VertexCode(G4String & aCode);

    void SetCode(G4String & aCode);
    G4String GetCode();
  private:

    G4String theCode;
};

inline void G4VertexCode::SetCode(G4String & aCode)
{
  theCode = aCode;
}

inline G4String G4VertexCode::GetCode()
{
  return theCode;
}

#endif
