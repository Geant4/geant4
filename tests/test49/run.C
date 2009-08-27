// ROOT Macro that launches the Geant4 model testing extensions
{	
	gSystem->Load("libGui.so");
	gSystem->Load("libG4ModelTester.so");

	G4TSimHelper::LoadLibraries();

	gGUIHelper->ShowMenu();
	//gGUIHelper->ShowPublicationGUI();
	//gGUIHelper->ShowSimulationGUI();
	//gGUIHelper->ShowAnalysisGUI();

	//gSimulationTool->Run("p_90_Al27_data.root", 450, "preco", 25, false);
	//gSimulationTool->Run("p_29_Al27_data.root", 740, "preco", 25, false);
	//gSimulationTool->Run("p_29_Au197_data.root", 1670, "preco", 25, false);

	
	//gAnalysisTool->Run("p_90_Al27_data.root",1);
	//gAnalysisTool->Run("p_29_Al27_data.root");
	//gAnalysisTool->Run("p_29_Au197_data.root");


}
