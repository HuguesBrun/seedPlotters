
TChain *chain = new TChain("eventsTree");

TH1F *giveEfficiency(TString nomPlot, TString variable, float *theBins, int nbBins){
    
    TH1F *denom = new TH1F("denom", "", nbBins, theBins);
    TH1F *num = new TH1F("num", "", nbBins, theBins);
    TString baseCut = "abs(T_Gen_Muon_PDGid)==13&&T_Gen_Muon_Pt>7&&abs(T_Gen_Muon_Eta)<2.4";
    chain->Draw(variable+">>denom",baseCut);
    chain->Draw(variable+">>num",baseCut+"&&T_Gen_Muon_L2crudeMaching==1");
    
    TH1F *efficiency = (TH1F*) denom->Clone(nomPlot);
    efficiency->Sumw2();
    efficiency->Divide(num, denom, 1,1);

    delete denom;
    delete num;
    return efficiency;
}



createEfficiencyPlots(TString dataset){
    chain->Add("../trees/tree_"+dataset+".root");
    
    TFile *myFile = new TFile("histoSingleMu_"+dataset+".root","RECREATE");

    
    float etaBins[11] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.1, 2.4};
    int nbEtaBins = 10;
    TH1F* histoEffEta = giveEfficiency("effVsEta", "abs(T_Gen_Muon_Eta)", etaBins, nbEtaBins);
    histoEffEta->Draw("E1");
    histoEffEta->Write();
    
    myFile->Close();
    delete myFile;
}
