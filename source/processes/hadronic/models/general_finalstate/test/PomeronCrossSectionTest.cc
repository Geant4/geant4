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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4PomeronCrossSection.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "globals.hh"
#include "G4Proton.hh"

main()
{
    G4Material *theCu = new G4Material("Copper", 8.96*g/cm3, 1);
    G4Element *elCu = new G4Element("Copper", "Cu",29.,63.55*g/mole);
    theCu->AddElement( elCu, 1 );
  double r0 = 5.0*fermi;
  int nStep = 10000;
  double db = r0/double(nStep);
  double s;
  
  int counter = 0;
  G4PomeronCrossSection * theProb = NULL;
  while(counter<1)
  {
    if(theProb==NULL) delete theProb;
    if(counter==0) theProb = new G4PomeronCrossSection(G4Proton::ProtonDefinition());
    if(counter==1) theProb = new G4PomeronCrossSection(G4Neutron::NeutronDefinition());
    if(counter==2) theProb = new G4PomeronCrossSection(G4PionPlus::PionPlusDefinition());
    if(counter==3) theProb = new G4PomeronCrossSection(G4PionMinus::PionMinusDefinition());
    if(counter==4) theProb = new G4PomeronCrossSection(G4KaonMinus::KaonMinusDefinition());
    if(counter==5) theProb = new G4PomeronCrossSection(G4KaonPlus::KaonPlusDefinition());
    theProb->Pomeron_S(2.7*GeV*GeV);
    theProb->Pomeron_Gamma((2.6+3.96)/GeV/GeV);
    theProb->Pomeron_C(1.4);
    theProb->Pomeron_Rsquare(3.56/GeV/GeV);
    theProb->Pomeron_Alpha(0.9808);
    theProb->Pomeron_Alphaprime(0.25/GeV/GeV);
    theProb->Pomeron_Gamma_Hard(0.0002/GeV/GeV);
    theProb->Pomeron_Alpha_Hard(1.47);
    double sqrtS = 5.*GeV;
    double errors;
    s = sqrtS*sqrtS;
    while(s<60*GeV*60*GeV)
    {
      s = sqrtS*sqrtS;
      double crossSection = 0;
      double inCrossSection = 0;
      double b = db/2.;
      for(int i=0; i<nStep; i++)
      {
    //    crossSection += 2.*3.14159265*b*db;
	crossSection += 2.*3.14159265*theProb->GetTotalProbability(s, b*b)*b*db;
	errors=theProb->GetTotalProbability(s,b*b);
	if(errors>1&&counter==2) G4cerr << s<<" "<<b<<" "<<errors<<G4endl;
        inCrossSection += 2.*3.14159265*theProb->GetInelasticProbability(s, b*b)*b*db;
	b+=db;
      }
      cout <<counter<<" "<< crossSection/millibarn<< " " <<sqrtS<<" ";
      cout << crossSection/millibarn-inCrossSection/millibarn<< " "<<G4endl;
  //    cout << "crossSection = "<< crossSection/millibarn<< " " <<sqrtS<< G4endl;
  //    cout << "inCrossSection = "<< inCrossSection/millibarn<< " "<<sqrtS<< G4endl;
  //    cout << 3.14159264*r0*r0/millibarn<< G4endl;
      sqrtS += 1.*GeV;
    }
    counter++;
  }
  delete theProb;
}
