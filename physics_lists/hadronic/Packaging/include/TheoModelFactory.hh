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
#ifndef TheoModelFactory_hh
#define TheoModelFactory_hh

template<class C, class S, class F>
class TheoModelFactory
{
  public:
  static G4TheoFSGenerator * New()
  {
    G4TheoFSGenerator * result = new G4TheoFSGenerator;
    G4VIntraNuclearTransportModel  * theCascade = new C;
    result->SetTransport(theCascade);

    G4VLongitudinalStringDecay * theFragmentation = new F;
    G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
    G4VPartonStringModel * theStringModel = new S;
    theStringModel->SetFragmentationModel(theStringDecay);
    result->SetHighEnergyGenerator(theStringModel);
    
    return result;
  }

  private:

};
// 2002 by J.P. Wellisch

#endif
