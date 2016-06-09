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
// $Id: G4FTFCrossSection.cc,v 1.2 2007/04/24 10:37:10 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

#include "G4FTFCrossSection.hh"

G4FTFCrossSection::G4FTFCrossSection()
{;}


G4FTFCrossSection::~G4FTFCrossSection()
{;}
//**********************************************************************************************

G4FTFCrossSection::G4FTFCrossSection(const G4ParticleDefinition * particle, G4double s) 
    {
    G4int PDGcode = particle->GetPDGEncoding();
    G4int absPDGcode = std::abs(PDGcode);
    G4double Elab = (s - 2*0.88*GeV*GeV)/(2*0.939*GeV)/GeV;
    G4double Plab = std::sqrt(Elab * Elab - 0.88);

    G4double LogPlab = std::log( Plab );
    G4double sqrLogPlab = LogPlab * LogPlab;

//G4cout<<"G4FTFCrossSection Plab "<<Plab<<G4endl;

    G4int NumberOfTargetProtons  = 1; //aNucleus.GetZ();             // ??????????????????????
    G4int NumberOfTargetNeutrons = 1; //aNucleus.GetN();
    G4int NumberOfTargetNucleons = NumberOfTargetProtons + NumberOfTargetNeutrons;

    G4double Xtotal, Xelastic;

    if( absPDGcode > 1000 )                        //------Projectile is baryon --------
      {        
       G4double XtotPP = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
       G4double XtotPN = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
       G4double XelPN  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotPP + 
                           NumberOfTargetNeutrons * XtotPN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelPP  + 
                           NumberOfTargetNeutrons * XelPN   ) / NumberOfTargetNucleons;
      }
    else if( PDGcode ==  211 )                     //------Projectile is PionPlus -------
      {
       G4double XtotPiP = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
       G4double XtotPiN = 33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
           
       G4double XelPiP  =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;
       G4double XelPiN  = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
      }
    else if( PDGcode == -211 )                     //------Projectile is PionMinus -------
      {
       G4double XtotPiP = 33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
       G4double XtotPiN = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
           
       G4double XelPiP  = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;
       G4double XelPiN  =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons;
      }

    else if( PDGcode ==  111 )                     //------Projectile is PionZero  -------
      {
       G4double XtotPiP =(16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab +   //Pi+
                          33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab)/2; //Pi-

       G4double XtotPiN =(33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab +   //Pi+
                          16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab)/2; //Pi-
           
       G4double XelPiP  =( 0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab +    //Pi+
                           1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-
       G4double XelPiN  =( 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab +   //Pi+
                           0.0  + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
      }
    else if( PDGcode == 321 )                      //------Projectile is KaonPlus -------
      {
       G4double XtotKP = 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab;
       G4double XtotKN = 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab;

       G4double XelKP  =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;
       G4double XelKN  =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
      }
    else if( PDGcode ==-321 )                      //------Projectile is KaonMinus ------
      {
       G4double XtotKP = 32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab - 5.6 *LogPlab;
       G4double XtotKN = 25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab - 2.9 *LogPlab;

       G4double XelKP  =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;
       G4double XelKN  =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
      }
    else if( PDGcode == 311 )                      //------Projectile is KaonZero ------
      {
       G4double XtotKP =( 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab +   //K+
                          32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab - 5.6 *LogPlab)/2; //K-
       G4double XtotKN =( 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab +   //K+
                          25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab - 2.9 *LogPlab)/2; //K-

       G4double XelKP  =(  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab +   //K+
                           7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab)/2; //K-
       G4double XelKN  =(  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab +   //K+
                           5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab)/2; //K-
       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
      }
    else                                          //------Projectile is undefined, Nucleon assumed
      {
       G4double XtotPP = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
       G4double XtotPN = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
       G4double XelPN  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotPP + 
                           NumberOfTargetNeutrons * XtotPN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelPP  + 
                           NumberOfTargetNeutrons * XelPN   ) / NumberOfTargetNucleons;
      };

      SetTotalCrossSection(Xtotal);
      SetElastisCrossSection(Xelastic);
      SetInelasticCrossSection(Xtotal-Xelastic);

//G4cout<<"G4FTFCrossSection Xt Xel "<<Xtotal<<" "<<Xelastic<<G4endl;

//-----------------------------------------------------------------------------------  
      SetSlope( Xtotal*Xtotal/16./pi/Xelastic/0.3894 ); // Slope parameter of elastic scattering
                                                        //      (GeV/c)^(-2))
                                                        // Gaussian parametrization of
                                                        // elastic scattering amplitude assumed

//G4cout<<"G4FTFCrossSection Slope "<<GetSlope()<<G4endl;

//-----------------------------------------------------------------------------------
      SetGamma0( GetSlope()*Xtotal/10./2./pi );
//-----------------------------------------------------------------------------------

//G4cout<<"G4FTFCrossSection Out"<<G4endl;

    } 
//**********************************************************************************************
