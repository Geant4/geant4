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
import hep.aida.*;
import java.util.Random;
import java.io.IOException;

/*
   Small plotter example
   compile with (provided aida is in your CLASSPATH variable)
      % javac plotResults.java
      % java -cp .:$CLASSPATH plotResult

*/


public class plotResults 
{
   public static void main(String[] argv) throws IOException
      {
	  String fileName="Pol01.aida";
	  if (argv.length>0) {
	      fileName=argv[0];
	      System.out.println("using file name "+fileName);
	  }

      IAnalysisFactory af = IAnalysisFactory.create();
      ITree tree = af.createTreeFactory().create(fileName,"xml");
      
      IHistogram1D h1 = (IHistogram1D) tree.find("1");
      IHistogram1D h2 = (IHistogram1D) tree.find("2");
      IHistogram1D h3 = (IHistogram1D) tree.find("3");
      IHistogram1D h4 = (IHistogram1D) tree.find("4");

      IPlotterFactory pf  = af.createPlotterFactory();
      IPlotter plotterPhoton = pf.create("Photon");
      plotterPhoton.createRegions(2,2);
      plotterPhoton.region(0).plot(h1);
      plotterPhoton.region(1).plot(h2);
      plotterPhoton.region(2).plot(h3);
      plotterPhoton.region(3).plot(h4);
      plotterPhoton.show();

      IHistogram1D h5 = (IHistogram1D) tree.find("5");
      IHistogram1D h6 = (IHistogram1D) tree.find("6");
      IHistogram1D h7 = (IHistogram1D) tree.find("7");
      IHistogram1D h8 = (IHistogram1D) tree.find("8");

      IPlotter plotterElectron = pf.create("Electron");
      plotterElectron.createRegions(2,2);
      plotterElectron.region(0).plot(h5);
      plotterElectron.region(1).plot(h6);
      plotterElectron.region(2).plot(h7);
      plotterElectron.region(3).plot(h8);
      plotterElectron.show();

      IHistogram1D h9 = (IHistogram1D) tree.find("9");
      IHistogram1D h10 = (IHistogram1D) tree.find("10");
      IHistogram1D h11 = (IHistogram1D) tree.find("11");
      IHistogram1D h12 = (IHistogram1D) tree.find("12");

      IPlotter plotterPositron = pf.create("Positron");
      plotterPositron.createRegions(2,2);
      plotterPositron.region(0).plot(h9);
      plotterPositron.region(1).plot(h10);
      plotterPositron.region(2).plot(h11);
      plotterPositron.region(3).plot(h12);
      plotterPositron.show();
   }
}
