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
// $Id: G4HadronValues.cc,v 1.16 2006/06/29 20:09:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

//
//  G4HadronValues class
//
//
//  Kinematic and dynamic values 
//  N.  Starkov 2003.
//
//  Modifications:
//  14.11.05 Use PDG code instead of static particle pointers (N.Starkov)
//  23.11.05 cleanup (V.Ivanchenko)
//

#include "globals.hh"
#include "G4HadronValues.hh"

G4HadronValues::G4HadronValues()
{}

G4HadronValues::~G4HadronValues()
{}

void 
G4HadronValues::GetHadronValues(const G4DynamicParticle* aHadron)
  {
       G4int iHadron(-1), iHadrCode;
       iHadrCode = aHadron->GetDefinition()->GetPDGEncoding();

//  G4cout<<" Code "<<iHadrCode<<G4endl;

       if(  iHadrCode == 2212 ||
            iHadrCode == 2112 ||
            iHadrCode == 3122 ||
            iHadrCode == 3222 ||
            iHadrCode == 3112 ||
            iHadrCode == 3212 ||
            iHadrCode == 3312 ||
            iHadrCode == 3322 ||
            iHadrCode == 3334 )   iHadron = 0;

       else if(
            iHadrCode == -2212 ||
            iHadrCode == -2112 ||
            iHadrCode == -3122 ||
            iHadrCode == -3222 ||
            iHadrCode == -3112 ||
            iHadrCode == -3212 ||
            iHadrCode == -3312 ||
            iHadrCode == -3322 ||
            iHadrCode == -3334 )   iHadron = 1;

       else if(  iHadrCode ==  211)     iHadron = 2;
       else if(  iHadrCode == -211)     iHadron = 3;
       else if(  iHadrCode ==  321)     iHadron = 4;
       else if(  iHadrCode == -321)     iHadron = 5;

       else {   
         G4cout << "G4HadronValues::GetHadronValues iHadrCode= " 
		<< iHadrCode
		<< "  " << aHadron->GetDefinition()->GetParticleName()
		<< G4endl;
         G4Exception(" There is not method for this hadron ");
       }

       G4double mHadr      = aHadron->GetMass()/1000.;         // In GeV
       G4double HadrEnergy = aHadron->GetTotalEnergy()/1000.;  // In GeV
       G4double HadrMoment = aHadron->GetTotalMomentum()/1000.;  // In GeV
       G4double sHadr      = 2*HadrEnergy*0.938+0.938*0.938+mHadr*mHadr;
       G4double sqrS       = std::sqrt(sHadr);
       G4double Ecm        = (sHadr-mHadr*mHadr+0.938*.938)/2/sqrS;
                MomentumCM = std::sqrt(Ecm*Ecm-0.938*0.938);

   if(HadrEnergy-mHadr<1.0) 
    {
     G4cout<<HadrEnergy<<G4endl;
     G4Exception(" The hadron Energy is very low for this method!");
    }

        switch (iHadron)
        {

         case 0:                  //  proton

        G4double Delta;

          Delta=1;

              if(HadrEnergy<40)
                  Delta = 0.916+0.0021*HadrEnergy;
              HadrTot   = 5.2+5.2*std::log(HadrEnergy)
                          +51*std::pow(HadrEnergy,-0.35);            //  mb
              HadrSlope = 6.44+0.88*std::log(sHadr)-1;               //  GeV-2 
              HadrReIm  = 0.13*std::log(sHadr/350)*std::pow(sHadr,-0.18);
              DDSect2   = 11;                                    //mb*GeV-2
              DDSect3   = 3;                                     //mb*GeV-2

//    if(HadrEnergy>1000) HadrReIm=0.15; 

        if( iHadrCode == 3122 || iHadrCode == 3222 ||
            iHadrCode == 3112 || iHadrCode == 3212 )
            {
              HadrTot   *=0.80;
              HadrSlope *=0.85;
            }

        if( iHadrCode == 3312 || iHadrCode == 3322 )
            {
              HadrTot   *=0.70;
              HadrSlope *=0.75;
            }           

         if( iHadrCode == 3334)
            {
              HadrTot   *=0.60;
              HadrSlope *=0.65;
            }

             break;

         case 1:              //   antiproton

              sqrS      = std::sqrt(sHadr);
              HadrTot   = 5.2+5.2*std::log(HadrEnergy)
                          +123.2*std::pow(HadrEnergy,-0.5);           //  mb
              HadrSlope = 8.32+0.57*std::log(sHadr); //GeV-2 
           if(HadrEnergy<1000)
              HadrReIm  =0.06*(sqrS-2.236)*(sqrS-14.14)*std::pow(sHadr,-0.8);
           else
              HadrReIm  = 0.6*std::log(sHadr/350)*std::pow(sHadr,-0.25);

              DDSect2   = 11;                                     //mb*GeV-2
              DDSect3   = 3;                                      //mb*GeV-2

//    if(HadrEnergy>1000) HadrReIm=0.15; 

        if( iHadrCode == -3122 || iHadrCode == -3222 ||
            iHadrCode == -3112 || iHadrCode == -3212 )
            {
              HadrTot   *=0.75;
              HadrSlope *=0.85;
            }

        if( iHadrCode == -3312 || iHadrCode == -3322 )
            {
              HadrTot   *=0.65;
              HadrSlope *=0.75;
            }

         if( iHadrCode == -3334)
            {
              HadrTot   *=0.55;
              HadrSlope *=0.65;
            }

           break;

         case 2:             //   pi plus

            if(HadrMoment>2.0)
              HadrTot    = 10.6+2.*std::log(HadrEnergy)+
                              25*std::pow(HadrEnergy,-0.43);           // mb
            else HadrTot = 40-50*(HadrMoment-1.5)*(HadrMoment-1.7);
              HadrSlope = 7.28+0.245*std::log(sHadr);                  //GeV-2 
              HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);
              DDSect2   = 4.6;                                    //mb*GeV-2
              DDSect3   = 1.33;                                   //mb*GeV-2
            break;

         case 3:             //   pi minus

              HadrTot   = 10.6+2*std::log(HadrEnergy)+
                          30*std::pow(HadrEnergy,-0.43);             // mb 

            if(HadrMoment<1.399)
              HadrTot = HadrTot+21.0/0.4*(1.4-HadrMoment);

              HadrSlope = 7.28+0.245*std::log(sHadr);               // GeV-2 
              HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);
              DDSect2   = 4.6;                                 //mb*GeV-2
              DDSect3   = 1.33;                                //mb*GeV-2
            break;

         case 4:            //  K plus

              HadrTot   = 10.6+1.8*std::log(HadrEnergy)+        
                               9.0*std::pow(HadrEnergy,-0.55);     // mb 
         if(HadrEnergy>100) HadrSlope = 15.0;
         else
//              HadrSlope = 5.28+1.76*std::log(sHadr)-
              HadrSlope = 1.0+1.76*std::log(sHadr)-
                              2.84*std::pow(sHadr,-0.5);           // GeV-2
              HadrReIm  = 0.4*(sHadr-20)*(sHadr-150)*std::pow(sHadr+50,-2.1);
              DDSect2   = 3.5;                                //mb*GeV-2
              DDSect3   = 1.03;                               //mb*GeV-2
            break;

         case 5:              //   K minus

              HadrTot   = 10+1.8*std::log(HadrEnergy)
                               +25*std::pow(HadrEnergy,-0.5);   // mb 
              HadrSlope = 6.98+0.127*std::log(sHadr);           // GeV-2 
//         if(HadrEnergy<8) HadrReIm = 0.7;
//         else
              HadrReIm  = 0.4*(sHadr-20)*(sHadr-20)*std::pow(sHadr+50,-2.1);
              DDSect2   = 3.5;                             //mb*GeV-2
              DDSect3   = 1.03;                            //mb*GeV-2
            break;
      }     
  }

/* end of file */
