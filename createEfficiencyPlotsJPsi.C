
TChain *chain = new TChain("eventsTree");

TH1F *giveEfficiency(TString nomPlot, TString variable, float *theBins, int nbBins){
    
    TH1F *denom = new TH1F("denom", "", nbBins, theBins);
    TH1F *num = new TH1F("num", "", nbBins, theBins);
    TString baseCut = "mu_pt>2&&abs(mu_eta)<2.4";
    chain->Draw(variable+">>denom",baseCut);
    chain->Draw(variable+">>num",baseCut+"&&foundL2==1");
    
    TH1F *efficiency = (TH1F*) denom->Clone(nomPlot);
    efficiency->Sumw2();
    efficiency->Divide(num, denom, 1,1);

    delete denom;
    delete num;
    return efficiency;
}



createEfficiencyPlotsJPsi(TString dataset){
    chain->Add("minitree_RelValJpsiMM_"+dataset+".root");
    
    TFile *myFile = new TFile("histoJPSI_"+dataset+".root","RECREATE");
    
    float etaBins[11] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.1, 2.4};
    int nbEtaBins = 10;
   TH1F* histoEffEta = giveEfficiency("effVsEta", "abs(mu_eta)", etaBins, nbEtaBins);
    histoEffEta->Write();
    
    float deltaRbins[8] = {0.2,0.5,0.6,0.7,0.8,1,1.5,5.0};
    int nbRbins = 7;
    TH1F* histoEffdR = giveEfficiency("effVsdR", "pair_dR", deltaRbins, nbRbins);
    histoEffdR->Write();

    
    float PtBinsbins[12] = {0, 4, 6, 8, 9, 10, 11, 12, 13, 14, 16, 20};
    int nbPtBins = 11;
    TH1F* histoEffPt = giveEfficiency("effVsPt", "pair_pt", PtBinsbins, nbPtBins);
    histoEffPt->Write();
   
    myFile->Close();
    delete myFile;
    
}
