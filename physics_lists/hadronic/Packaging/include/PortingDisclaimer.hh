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
#ifndef PortingDisclaimer_h
#define PortingDisclaimer_h

class PortingDisclaimer
{
  protected:
  PortingDisclaimer()
  {
      G4Exception("PortingDisclaimer", "007", FatalException,
                   "The disclaimer is not to be instanciated" );
  }
  private:
  struct LocalDisclaimer
  {
    LocalDisclaimer()
    {
      G4cout << G4endl;
      G4cout << "========================== porting disclaimer =================================="<<G4endl;
      G4cout << "You are using a technical port of the physics lists provided"<<G4endl
             << "for hadronic use-cases, as validated for geant4 5.2, with a "<<G4endl
	     << "later version of geant4."<<G4endl
	     << "These physics lists are made available with that geant4 districution, so users"<<G4endl
	     << "wishing to use the latest-greatest geant4 code keep having a starting point"<<G4endl
	     << "for their work with hadronics."<<G4endl<<G4endl
	     << "Please note that it is a purely technical port, without physics validation performed. "<<G4endl
	     << "Hence these physics lists do not yet come with the usual guarantees attached as to "<<G4endl
	     << "their physics performance."<< G4endl
	     << "For validated 'educated guess' physics lists for hadronic usecases, "<<G4endl
	     << "as well as updates, please see the related WWW pages."<< G4endl << G4endl
	     << "For information as to which physics list is optimized for your study,"<<G4endl
	     << "please also consult the GHAD WWW pages."<<G4endl    
	     << "Geant4 hadronics is a component system, and here the physics lists main "<<G4endl
	     << "task is to configure these components in a suitable way for a concrete useage."<<G4endl
	     << "Not all physics lists are suitable for all usecases."<<G4endl;        
      G4cout << "========================== porting disclaimer =================================="<<G4endl;
      G4cout << G4endl;
    }
  };
  public:
  static LocalDisclaimer * PrintLocalDisclaimer() 
  {
    static LocalDisclaimer aD;
    return &aD;
  }
};

#endif
