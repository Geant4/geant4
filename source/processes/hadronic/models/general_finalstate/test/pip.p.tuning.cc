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
  double db = r0/double(nStep); // (24, 51, 480)
  double s;
  const double ss[] = {7.*GeV, 10.*GeV, 30.*GeV};
  const double xtot[] = {23.5*millibarn, 23.0*millibarn, 24.8*millibarn};
  const double xel[] = {3.5*millibarn, 3.29*millibarn, 3.15*millibarn};
  
  int counter = 0;
  G4PomeronCrossSection * theProb = NULL;
  theProb = new G4PomeronCrossSection(G4Proton::ProtonDefinition());
  G4bool done = false;
  double oldchi2 = 0;
  
  double dPomeron_S = 0.1*GeV*GeV;
  double epsPomeron_S = 0.*GeV*GeV;
  
  double dPomeron_Gamma = 0.02/GeV/GeV;
  double epsPomeron_Gamma = 1.54/GeV/GeV;
  
  double dPomeron_C = 0.1;
  double epsPomeron_C = 0.;
  
  double dPomeron_Rsquare = 0.1/GeV/GeV;
  double epsPomeron_Rsquare = 0./GeV/GeV;
  
  double dPomeron_Alpha = 0.025;
  double epsPomeron_Alpha = -0.1;
  
  double dPomeron_Alphaprime = 0.1/GeV/GeV;
  double epsPomeron_Alphaprime = 0./GeV/GeV;
  
  double dPomeron_Gamma_Hard = 0.1/GeV/GeV;
  double epsPomeron_Gamma_Hard = 0./GeV/GeV;
  
  double dPomeron_Alpha_Hard = 0.1;
  double epsPomeron_Alpha_Hard = 0;
  
  while(!done)
  {
    double chi2 = 0;
    double dchi2;
    counter = 0;
    double crossSection = 0;
    double inCrossSection = 0;
    while(counter<3)
    {
      theProb->Pomeron_S(1.5*GeV*GeV + epsPomeron_S);
      theProb->Pomeron_Gamma(2.17/GeV/GeV + epsPomeron_Gamma);
      theProb->Pomeron_C(1.6 + epsPomeron_C);
      theProb->Pomeron_Rsquare(2.36/GeV/GeV + epsPomeron_Rsquare);
      theProb->Pomeron_Alpha(1.0808 + epsPomeron_Alpha);
      theProb->Pomeron_Alphaprime(0.25/GeV/GeV + epsPomeron_Alphaprime);
      theProb->Pomeron_Gamma_Hard(0.0002/GeV/GeV + epsPomeron_Gamma_Hard);
      theProb->Pomeron_Alpha_Hard(1.47 + epsPomeron_Alpha_Hard);
      s = ss[counter]*ss[counter];
      crossSection = 0;
      inCrossSection = 0;
      double b = db/2.;
      for(int i=0; i<nStep; i++)
      {
	crossSection += 2.*3.14159265*theProb->GetTotalProbability(s, b*b)*b*db;
	inCrossSection += 2.*3.14159265*theProb->GetInelasticProbability(s, b*b)*b*db;
	b+=db;
      }
      dchi2 = crossSection-xtot[counter];
      chi2 += dchi2*dchi2;
      dchi2 = crossSection-inCrossSection-xel[counter];
      chi2 += dchi2*dchi2;
      counter++;
    cout << "Xsec = "<<crossSection<<", Xel = "<<
             crossSection-inCrossSection<<G4endl;
    }
//    epsPomeron_S += dPomeron_S; // unknown
//    epsPomeron_Gamma += dPomeron_Gamma; // absolute norm
//    epsPomeron_C += dPomeron_C; // abs norm of elastic
//    epsPomeron_Rsquare += dPomeron_Rsquare; // unknown
//    epsPomeron_Alpha -= dPomeron_Alpha; // even more interesting - elastic slope
    epsPomeron_Alphaprime -= dPomeron_Alphaprime; // interesting for elastic X-sec
//    epsPomeron_Gamma_Hard += dPomeron_Gamma_Hard;
//    epsPomeron_Alpha_Hard += dPomeron_Alpha_Hard;
    cout << "Xsec = "<<crossSection<<", Xel = "<<
             crossSection-inCrossSection<<", chi2 = "<<chi2<<G4endl;
    cout << "for "<<epsPomeron_Gamma<< "=================================================="<<G4endl;
  }
  delete theProb;
}
