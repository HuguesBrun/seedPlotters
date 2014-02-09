

drawThePlot(TString  nomFile, TString nomHisto, TString titleName){
    
    TFile *myFileOld = new TFile("histoJPSI_oldGeom.root");
    TFile *myFileNew = new TFile("histoJPSI_PostLS1.root");
    
    TH1F *h_old  = (TH1F*) myFileOld->Get(nomHisto);
    TH1F *h_New  = (TH1F*) myFileNew->Get(nomHisto);
    
    
    TCanvas *c0 = new TCanvas("c0","coucou",600,600);
    c0->SetFillColor(0);
    h_old->SetFillColor(kGreen-7);
    h_old->GetXaxis()->SetTitle(titleName);;
    h_old->GetYaxis()->SetTitle("efficiency");
    h_old->SetMaximum(1.1);
    h_old->SetMinimum(0.5);
    h_old->Draw("E3");
    h_New->SetFillColor(kMagenta-7);
    h_New->SetLineWidth(3);
    h_New->SetLineColor(kMagenta-7);
    h_New->Draw("E1:same");
    
    TLegend *t = new TLegend(0.57,0.14,0.85,0.31);
    t->SetFillColor(0);
    t->SetLineColor(0);
    t->AddEntry(h_old,"run1 geometry","F");
    t->AddEntry(h_New,"postLS1 geometry","L");
    t->Draw();
    
    c0->Print("plotsJPsi/"+nomFile+".png");
    
    delete c0;
    
}

drawHistoJPsi(){
    gStyle->SetOptStat(0);
    drawThePlot("plot_effVsEta", "effVsEta", "gen muon |#eta|");
    drawThePlot("plot_effVsdR", "effVsdR", "gen #Delta R from closest neighbour");
    drawThePlot("plot_effVsPt", "effVsPt", "gen P_{T} of J/#psi");

}
