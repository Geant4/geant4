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
//
// $Id: G4QPomeron.cc,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QPomeron ----------------
//      by Mikhail Kossov Oct, 2006
//      class for a Pomeron used by Parton String Models
//   For comparison mirror member functions are taken from G4 class:
//   G4PomeronCrossSection
// -----------------------------------------------------------------

#include "G4QPomeron.hh"

G4QPomeron::G4QPomeron(G4int PDG)
{
	 pomeron_Alpha=		1.0808;
	 pomeron_Alphaprime =	0.25/GeV/GeV;     
  pomeron_Gamma_Hard = 0.0002/GeV/GeV;
  pomeron_Alpha_Hard = 1.47;
  G4int aP = std::abs(PDG);
  if     (aP==2212||aP==2112)                   InitForNucleon();
  else if(PDG==111||aP==211)                    InitForPion();
  else if(PDG==130||PDG==310||aP==311||aP==321) InitForKaon();
  else if(PDG==22)                              InitForGamma();
  else if(aP==3122||aP==3222||aP==3212||aP==3112||aP==3322||aP==3312||aP==3334)
                                                InitForNucleon();
  else
		{
    G4cout<<"-Warning-G4QPomeron is initialised for PDGCode="<<PDG<<" as a pion"<<G4endl;
    InitForPion();
  }
}

G4double G4QPomeron::GetCutPomeronProbability(const G4double s,const G4double imp2,
                                              const G4int nPom)
{
	 G4double f=1; if(nPom>1) for(G4int i=2; i<= nPom; i++) f*=i; // Calculate factorial
  G4double e=Eikonal(s,imp2); e+=e;                            // Doubled Eikonal
 	return std::exp(-e)*std::pow(e,nPom)/pomeron_C/f;
}

void G4QPomeron::InitForNucleon()
{
  // pomeron_S=		3.0*GeV*GeV;  
	 pomeron_S=		2.7*GeV*GeV;
  //	pomeron_Gamma=		2.16/GeV/GeV;      // ? M.K.
  //	pomeron_Gamma=		3.96/GeV/GeV;      // ? M.K.
	 pomeron_Gamma   =	(2.6+3.96)/GeV/GeV; // ? M.K.
	 pomeron_C       =	1.4;
	 pomeron_Rsquare =	3.56/GeV/GeV;
	 /////pomeron_Alpha=		0.9808;          // This is very strange M.K. (was in GHAD QGSM)
  // If pomeron_Gamma_Hard!=0 to fit total pp XS use pomeron_Gamma_Soft=2.35/GeV/GeV ? M.K.
}

void G4QPomeron::InitForPion()
{
	 pomeron_S       =	1.5*GeV*GeV;
  //pomeron_Gamma =	1.46/GeV/GeV;  // ? M.K.
	 pomeron_Gamma   =	2.17/GeV/GeV;  // ? M.K.
	 pomeron_C       =	1.6;
	 pomeron_Rsquare =	2.36/GeV/GeV;
}

void G4QPomeron::InitForKaon()
{
	 pomeron_S       =	2.3*GeV*GeV;
  //pomeron_Gamma =	1.31/GeV/GeV;  // ? M.K.
	 pomeron_Gamma   = 1.92/GeV/GeV;  // ? M.K.
	 pomeron_C       =	1.8;
	 pomeron_Rsquare =	1.96/GeV/GeV;
}

void G4QPomeron::InitForGamma()
{
	 pomeron_S       =	1.7*GeV*GeV;
  //pomeron_Gamma =		1.42/GeV/GeV; // ? M.K.
	 pomeron_Gamma   =	2.07/GeV/GeV;  // ? M.K.
	 pomeron_C       =	1.7;
	 pomeron_Rsquare =	2.16/GeV/GeV;
}

G4double G4QPomeron::Expand(G4double z)
{

	G4double sum=1.;
	G4double current=1.;
	for(G4int j=2; j<21; j++)
	{
		 current *= -z*(j-1)/j/j;
		 sum+=current;
	}
	return sum;
}








    


