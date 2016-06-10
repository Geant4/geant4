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
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
//
// $Id: G4SmpNEngCf252.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"

G4double G4fissionEvent::G4SmpNEngCf252(G4int option) {

/*
  Description
    Sample energy spectrum for Cf-252. 
*/

/*
  Input
    option    - 0 Mannhart corrected Maxwellian spectrum
                1 Madland-Nix theoretical spectrum
                2 Froehner Watt spectrum
  Return
    energy of neutron emitted by spontaneous fission
*/

  G4double a,b;
  G4double gpar;
  G4double g2;
  G4double ferg;
  G4double r;

  r = fisslibrng();

/*
   Mannhart Corrected Spectrum
*/
  if(option == 0) {
    if(r == 0) return 0.000001;
    if (r > 0.0 && r <= 0.0005)
      return 0.00003 + 0.04992*(r/0.0005) - 0.59473*std::pow(r/0.0005,2) 
           + 5.44877*std::pow(r/0.0005,3) - 29.38086*std::pow(r/0.0005,4)
           + 97.14014*std::pow(r/0.0005,5) - 202.82112*std::pow(r/0.0005,6)
           + 268.2301*std::pow(r/0.0005,7) - 217.75316*std::pow(r/0.0005,8)
           + 98.96285*std::pow(r/0.0005,9) - 19.27077*std::pow(r/0.0005,10);

    if (r > 0.0005 && r <= 0.005)
      return 0.01118 + 0.06715*((r-.0005)/.0045)
           - 0.09236*std::pow((r-.0005)/.0045,2) + 0.26224*std::pow((r-.0005)/.0045,3)
           - 0.64784*std::pow((r-.0005)/.0045,4) + 1.16830*std::pow((r-.0005)/.0045,5)
           - 1.43858*std::pow((r-.0005)/.0045,6) + 1.13771*std::pow((r-.0005)/.0045,7)
           - 0.51839*std::pow((r-.0005)/.0045,8) + 0.10302*std::pow((r-.0005)/.0045,9);

    if (r > 0.005 && r <= 0.05)
      return 0.05244+0.32101*((r-.005)/.045)
           -  0.52574*std::pow((r-.005)/.045,2) +   2.80540*std::pow((r-.005)/.045,3)
           - 14.88036*std::pow((r-.005)/.045,4) +  55.46869*std::pow((r-.005)/.045,5)
           -133.64517*std::pow((r-.005)/.045,6) + 202.88434*std::pow((r-.005)/.045,7)
           -186.86758*std::pow((r-.005)/.045,8) +  95.19530*std::pow((r-.005)/.045,9)
           - 20.55275*std::pow((r-.005)/.045,10);

    if(r > 0.05 && r <= 0.25) return 0.25585+0.75532*((r-.05)/.2)-0.73676*std::pow((r-.05)/.2,2)+3.65653*std::pow((r-.05)/.2,3)-13.80528*std::pow((r-.05)/.2,4)+33.35932*std::pow((r-.05)/.2,5)-50.0410*std::pow((r-.05)/.2,6)+45.13793*std::pow((r-.05)/.2,7)-22.4072*std::pow((r-.05)/.2,8)+4.70141*std::pow((r-.05)/.2,9);

    if(r > 0.25 && r <= 0.50) return 0.87609+0.74687*((r-.25)/.25)+0.02849*std::pow((r-.25)/.25,2)+0.06145*std::pow((r-.25)/.25,3)-0.09589*std::pow((r-.25)/.25,4)+0.29798*std::pow((r-.25)/.25,5)-0.57707*std::pow((r-.25)/.25,6)+0.66181*std::pow((r-.25)/.25,7)-0.40720*std::pow((r-.25)/.25,8)+0.10370*std::pow((r-.25)/.25,9);

    if(r > 0.5 && r <= 0.75) return 1.69622+0.93896*((r-.5)/.25)+0.16428*std::pow((r-.5)/.25,2)+0.21761*std::pow((r-.5)/.25,3)-0.96904*std::pow((r-.5)/.25,4)+3.34951*std::pow((r-.5)/.25,5)-6.35177*std::pow((r-.5)/.25,6)+6.90120*std::pow((r-.5)/.25,7)-3.98682*std::pow((r-.5)/.25,8)+0.95276*std::pow((r-.5)/.25,9);

    if(r > 0.75 && r <= 0.95) return 2.91217+1.52474*((r-.75)/.2)-4.99340*std::pow((r-.75)/.2,2)+58.72977*std::pow((r-.75)/.2,3)-313.30984*std::pow((r-.75)/.2,4)+946.0791*std::pow((r-.75)/.2,5)-1679.85559*std::pow((r-.75)/.2,6)+1740.83984*std::pow((r-.75)/.2,7)-973.51886*std::pow((r-.75)/.2,8)+227.06831*std::pow((r-.75)/.2,9);
     if(r > 0.95 && r <= 0.975) return 5.50137-0.99765*((r-.95)/.025)+27.57678*std::pow((r-.95)/.025,2)-218.47931*std::pow((r-.95)/.025,3)+1024.0426*std::pow((r-.95)/.025,4)-3005.86182*std::pow((r-.95)/.025,5)+5684.52295*std::pow((r-.95)/.025,6)-6919.36182*std::pow((r-.95)/.025,7)+5235.71777*std::pow((r-.95)/.025,8)-2240.06934*std::pow((r-.95)/.025,9)+413.9299*std::pow((r-.95)/.025,10);

     if(r > 0.975 && r <= 0.995) return 6.52172+1.21273*((r-.975)/.02)+0.69998*std::pow((r-.975)/.02,2)-1.78886*std::pow((r-.975)/.02,3)+11.57883*std::pow((r-.975)/.02,4)-39.41592*std::pow((r-.975)/.02,5)+88.32992*std::pow((r-.975)/.02,6)-127.68685*std::pow((r-.975)/.02,7)+115.97678*std::pow((r-.975)/.02,8)-60.09069*std::pow((r-.975)/.02,9)+13.66798*std::pow((r-.975)/.02,10);
     if(r > 0.995 && r <= 0.999) return 9.00502+1.31798*((r-.995)/.004)-1.17448*std::pow((r-.995)/.004,2)+20.15941*std::pow((r-.995)/.004,3)-114.27763*std::pow((r-.995)/.004,4)+370.04855*std::pow((r-.995)/.004,5)-701.888*std::pow((r-.995)/.004,6)+776.28204*std::pow((r-.995)/.004,7)-462.68823*std::pow((r-.995)/.004,8)+115.05296*std::pow((r-.995)/.004,9);
     if(r > 0.999 && r <= 0.9997) return 11.83792-1.8952*((r-.999)/.0007)+50.30901*std::pow((r-.999)/.0007,2)-239.56978*std::pow((r-.999)/.0007,3)+514.90747*std::pow((r-.999)/.0007,4)-508.73672*std::pow((r-.999)/.0007,5)+191.09637*std::pow((r-.999)/.0007,6);
     if(r > 0.9997) return 20.;
  }
/*
   Madland-Nix Spectrum
*/
  if(option == 1) {
     if(r <= 1.001065092e-03) return 1.946313876*std::pow(r,0.6667261950);
     else if(r > 1.001065092e-03 && r <= 1.001389105e-02) return 2.00504119*std::pow(r,0.6709990736);
     else if(r > 1.001389105e-02 && r <= 5.022359145e-02) return 2.107978578*std::pow(r,0.7077041191);
     else if(r > 5.022359145e-02 && r <= 1.000989427e-01) return 2.280517358*std::pow(r,0.7077041191);
     else if(r > 1.000989427e-01 && r <= 1.500872491e-01) return 2.444108408*std::pow(r,0.73764526215);
     else if(r > 1.500872491e-01 && r <= 2.002079974e-01) return 2.621855634*std::pow(r,0.7745779546);
     else if(r > 2.002079974e-01 && r <= 2.25221648e-01) return 2.753099265*std::pow(r,0.8044994010);
     else if(r > 2.25221648e-01 && r <= 2.501564538e-01) return 2.834010751*std::pow(r,0.8239187384);
     else if(r > 2.501564538e-01 && r <= 2.752546770e-01) return 2.911676280*std::pow(r,0.8434235719);
     else if(r > 2.752546770e-01 && r <= 3.000964724e-01) return 2.988430135*std::pow(r,0.8635883266);
     else if(r > 3.000964724e-01 && r <= 3.500470095e-01) return 3.099471293*std::pow(r,0.8942289512);
     else if(r > 3.500470095e-01 && r <= 4.001118970e-01) return 3.244686176*std::pow(r,0.9378302608);
     else if(r > 4.001118970e-01 && r <= 5.000461778e-01) return 3.543403932*std::pow(r,1.0411008510);
     else if(r > 5.000461778e-01 && r <= 5.501318506e-01) return 3.708358099*std::pow(r,1.1068317830);
     else if(r > 5.501318506e-01 && r <= 6.000655433e-01) return 3.889805304*std::pow(r,1.1868908770);
     else if(r > 6.000655433e-01 && r <= 6.500147305e-01) return 4.092497225*std::pow(r,1.2865658570);
     else if(r > 6.500147305e-01 && r <= 7.000271284e-01) return 4.322906068*std::pow(r,1.4140909190);
     else if(r > 7.000271284e-01 && r <= 7.501159110e-01) return 4.589909069*std::pow(r,1.5828217210);
     else if(r > 7.501159110e-01 && r <= 8.000662513e-01) return 4.906598744*std::pow(r,1.8162034790);
     else if(r > 8.000662513e-01 && r <= 8.500772033e-01) return 5.297053797*std::pow(r,2.1626825870);
     else if(r > 8.500772033e-01 && r <= 8.750123088e-01) return 5.650277904*std::pow(r,2.5517142900);
     else if(r > 8.750123088e-01 && r <= 9.000106866e-01) return 5.947741976*std::pow(r,2.9383159800);
     else if(r > 9.000106866e-01 && r <= 9.250286977e-01) return 6.317014169*std::pow(r,3.5155713570);
     else if(r > 9.250286977e-01 && r <= 9.350074655e-01) return 6.625757778*std::pow(r,4.1118364020);
     else if(r > 9.350074655e-01 && r <= 9.400070002e-01) return 6.784126941*std::pow(r,4.4594479870);
     else if(r > 9.400070002e-01 && r <= 9.500026229e-01) return 6.969180156*std::pow(r,4.9019105900);
     else if(r > 9.500026229e-01 && r <= 9.600065896e-01) return 7.254643542*std::pow(r,5.6894827520);
     else if(r > 9.600065896e-01 && r <= 9.700165577e-01) return 7.613500497*std::pow(r,6.8841593900);
     else if(r > 9.700165577e-01 && r <= 9.750157135e-01) return 7.944100103*std::pow(r,8.2544400860);
     else if(r > 9.750157135e-01 && r <= 9.800101585e-01) return 8.228439642*std::pow(r,9.6531190300);
     else if(r > 9.800101585e-01 && r <= 9.850018119e-01) return 8.586524083*std::pow(r,11.783756400);
     else if(r > 9.850018119e-01 && r <= 9.875072929e-01) return 8.917364901*std::pow(r,14.240137310);
     else if(r > 9.875072929e-01 && r <= 9.900006975e-01) return 9.202675761*std::pow(r,16.76089029);
     else if(r > 9.900006975e-01 && r <= 9.925048152e-01) return 9.562781386*std::pow(r,20.61962568);
     else if(r > 9.925048152e-01 && r <= 9.935030103e-01) return 9.867915664*std::pow(r,24.69147261);
     else if(r > 9.935030103e-01 && r <= 9.945000177e-01) return 10.08727342*std::pow(r,28.07701487);
     else if(r > 9.945000177e-01 && r <= 9.950025127e-01) return 10.27382614*std::pow(r,31.36001051);
     else if(r > 9.950025127e-01 && r <= 9.955029368e-01) return 10.41724243*std::pow(r,34.13127669);
     else if(r > 9.955029368e-01 && r <= 9.960005970e-01) return 10.57636221*std::pow(r,37.50088614);
     else if(r > 9.960005970e-01 && r <= 9.965016080e-01) return 10.75639015*std::pow(r,41.72354164);
     else if(r > 9.965016080e-01 && r <= 9.970001795e-01) return 10.96366661*std::pow(r,47.18729543);
     else if(r > 9.970001795e-01 && r <= 9.975004375e-01) return 11.20771170*std::pow(r,54.54899604);
     else if(r > 9.975004375e-01 && r <= 9.978504408e-01) return 11.45202216*std::pow(r,63.11906699);
     else if(r > 9.978504408e-01 && r <= 9.989524675e-01)
       return 2.72756636666e5-5.47258138432e5*r+2.74514044871e5*std::pow(r,2);
     else if(r > 9.989524675e-01 && r <= 9.994929298e-01)
       return 1.14946879661e6-2.30252188973e6*r+1.15306661788e6*std::pow(r,2);
     else if(r > 9.994929298e-01 && r <= 9.997558922e-01)
       return 4.90621526236e6-9.81982943883e6*r+4.91362868673e6*std::pow(r,2);
     else if(r > 9.997558922e-01 && r <= 9.998830120e-01)
       return 2.11365688795184e7-4.22884732250404e7*r+2.11519198434219e7*std::pow(r,2);
     else if(r > 9.998830120e-01 && r <= 9.999441620e-01)
       return 9.18987945911229e7-1.83829506875257e8*r+9.19307287711182e7*std::pow(r,2);
     else if(r > 9.999441620e-01 && r <= 9.999734440e-01)
       return 4.02781481130433e8-8.05629656768407e8*r+4.02848193115356e8*std::pow(r,2);
     else if(r > 9.999734440e-01 && r <= 9.999874120e-01)
       return 1.77804635135775e9-3.55623257045546e9*r+1.77818623756641e9*std::pow(r,2);
     else if(r > 9.999874120e-01 && r <= 9.999940510e-01)
       return 7.90099032702915e9-1.58022749659903e10*r+7.90128465842187e9*std::pow(r,2);
     else if(r > 9.999940510e-01 && r <= 9.999971960e-01)
       return 3.53223507413091e10-7.06453227162775e10*r+3.53229719954219e10*std::pow(r,2);
     else if(r > 9.999971960e-01 && r <= 9.999986820e-01)
       return 1.58786475903785e11-3.17574266841213e11*r+1.58787790958875e11*std::pow(r,2);
     else if(r > 9.999986820e-01 && r <= 9.999993820e-01)
       return 7.17433904438156e11-1.43487059972047e12*r+7.17436695304750e11*std::pow(r,2);
     else if(r > 9.999993820e-01 && r <= 9.999997110e-01)
       return 3.257374123945330e12-6.514754184993900e12*r+3.257380061072000e12*std::pow(r,2);
     else if(r > 9.999997110e-01 && r <= 9.999998650e-01)
       return 1.48641255466171e13-2.97282637539286e13*r+1.48641382073360e13*std::pow(r,2);
     else if(r > 9.999998650e-01 && r <= 9.999999370e-01)
       return 6.82056055248876e13-1.36411238119518e14*r+6.82056325946560e13*std::pow(r,2);
     else if(r > 9.999999370e-01 && r <= 1.000000000e00)
       return 3.14919363013517e14-6.29838784079090e14*r+3.14919421065600e14*std::pow(r,2);
  }
/*
   Frohner Watt Spectrum
*/
  if (option == 2) {
    a=1.175;
    b=1.040;

    do {
      gpar = std::sqrt(std::pow(1+0.125*a*b,2.)-1)+(1+0.125*a*b);
      g2=-std::log(fisslibrng());
      ferg=a*gpar*g2;
    } while (std::pow((1-gpar)*(1+g2)-std::log(fisslibrng()),2.) > b*ferg);
    return ferg;
  }

  //
  // Fall through
  //
   
  G4cout << " SmpNEngCf252: unrecognized option = " << option << G4endl;
  return -1.0;
}
