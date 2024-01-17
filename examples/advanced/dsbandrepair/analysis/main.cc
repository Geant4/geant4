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
/// \file main.cc
/// \brief Main program of the Analysis module

#include "AnalysisHandler.hh"
#include <iostream>
#include "ParametersParser.hh"
#include "SDDData.hh"
int main(int argc,char** argv)
{
    std::cout  <<"#####################################################################\n"
            <<"#                            dsbandrepair                           #\n"
            <<"#                Welcome to \"Analysis Module\" v.1.0                 #\n"
            <<"#####################################################################\n"
            <<"\n"
            <<"--------------------------> Start running <--------------------------"<<std::endl;
    ParametersParser *parParser = ParametersParser::Instance();
    if (argc > 1) {
        std::string macrofile = argv[1];
        parParser->LoadParameters(macrofile);
    }
    
    AnalysisHandler aAna;
    if (parParser->GetBpForDSB() > 0) aAna.SetBpForDSB(parParser->GetBpForDSB());
    aAna.GiveMeSBs();
    if (!parParser->WannaLoadDamagesFromSDD()) {
        aAna.CreateSDD(("SDDformat_"+parParser->GetOutputName()));
    }
    
    if (parParser->UseTLK()) aAna.ApplyDNAModel("TLK");
    if (parParser->UseLEMIV()) aAna.ApplyDNAModel("LEMIV");
    if (parParser->UseBelov()) aAna.ApplyDNAModel("BELOV");
    std::cout  <<"----------------------> Finish!!! Good bye :) <----------------------"<<std::endl;
    return 0;
}
