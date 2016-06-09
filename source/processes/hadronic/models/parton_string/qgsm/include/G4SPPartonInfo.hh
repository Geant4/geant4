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
#ifndef G4SPPartonInfo_h
#define G4SPPartonInfo_h

class G4SPPartonInfo
{
  public:
    G4SPPartonInfo(G4int diq, G4int q, G4double prob) 
    { diQuarkPDGCode = diq; quarkPDGCode = q; probability = prob; }
    G4int GetQuark() const {return quarkPDGCode;}
    G4int GetDiQuark() const {return diQuarkPDGCode;}
    G4double GetProbability() const {return probability;}      
    G4bool operator == (const G4SPPartonInfo & aInfo) const
    {return this == &aInfo;}
  private:      
    G4int quarkPDGCode;
    G4int diQuarkPDGCode;
    G4double probability;
};
    
#endif
