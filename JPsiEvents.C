TString nomFichier = "RelValJpsiMM_OldGeom";


TChain *chain = new TChain("eventsTree");

std::vector<float> *T_Gen_Muon_Pt;
std::vector<float> *T_Gen_Muon_Eta;
std::vector<float> *T_Gen_Muon_Phi;
std::vector<float> *T_Gen_Muon_Energy;
std::vector<float> *T_Gen_Muon_Px;
std::vector<float> *T_Gen_Muon_Py;
std::vector<float> *T_Gen_Muon_Pz;
std::vector<float> *T_Gen_Muon_tpPt;
std::vector<float> *T_Muon_Eta;
std::vector<float> *T_Muon_Phi;
std::vector<float> *T_Muon_Pt;
std::vector<bool> *T_Muon_IsGlobalMuon;
std::vector<bool> *T_Muon_IsTrackerMuon;
std::vector<int> *T_Gen_Muon_PDGid;
std::vector<int> *T_Gen_Muon_status;
std::vector<int> *T_Gen_Muon_MotherID;
std::vector<int> *T_Gen_Muon_FoundSTA;
std::vector<int> *T_Gen_Muon_FoundL2;
std::vector<int> *T_Gen_Muon_L2crudeMaching;


TFile *myOutFile = new TFile("minitree_"+nomFichier+".root","RECREATE");

TTree *mytree;
float mu_eta;
float mu_pt;
int   foundL2;
int   foundL2crude;
int   foundSTA;
bool   foundRECO;
bool   foundRECOglobal;
bool   foundRECOtracker;
float pair_mass;
float pair_dR;
float pair_pt;


float deltaR(float phi1, float phi2, float eta1, float eta2)
{
    float dphi=deltaPhi(phi1,phi2);
    float deta=fabs(eta1-eta2);
    float dr = sqrt(dphi*dphi+ deta*deta);
    return dr;
}

float deltaPhi(float phi1, float phi2)
{
    float dphi;
    if(phi1<0) phi1+=2*TMath::Pi();
    if(phi2<0) phi2+=2*TMath::Pi();
    dphi=fabs(phi1-phi2);
    if(dphi>2*TMath::Pi()) dphi-=2*TMath::Pi();
    if(dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
    return dphi;
}



JPsiEvents(){
    
    mytree = new TTree("eventsTree","");
    mytree->Branch("mu_eta", &mu_eta, "mu_eta/F");
    mytree->Branch("mu_pt", &mu_pt, "mu_pt/F");
    mytree->Branch("pair_mass", &pair_mass, "pair_mass/F");
    mytree->Branch("pair_pt", &pair_pt, "pair_pt/F");
    mytree->Branch("pair_dR", &pair_dR, "pair_dR/F");
    mytree->Branch("foundL2", &foundL2, "foundL2/I");
    mytree->Branch("foundL2crude", &foundL2crude, "foundL2crude/I");
    mytree->Branch("foundSTA", &foundSTA, "foundSTA/I");
    mytree->Branch("foundRECO", &foundRECO, "foundRECO/B");
    mytree->Branch("foundRECOglobal", &foundRECOglobal, "foundRECOglobal/B");
    mytree->Branch("foundRECOtracker", &foundRECOtracker, "foundRECOtracker/B");
    
    chain->Add("../trees/tree_"+nomFichier+".root");
    chain->SetBranchAddress("T_Gen_Muon_Eta",&T_Gen_Muon_Eta);
    chain->SetBranchAddress("T_Gen_Muon_Phi",&T_Gen_Muon_Phi);
    chain->SetBranchAddress("T_Gen_Muon_Pt",&T_Gen_Muon_Pt);
    chain->SetBranchAddress("T_Gen_Muon_Energy",&T_Gen_Muon_Energy);
    chain->SetBranchAddress("T_Gen_Muon_Px",&T_Gen_Muon_Px);
    chain->SetBranchAddress("T_Gen_Muon_Py",&T_Gen_Muon_Py);
    chain->SetBranchAddress("T_Gen_Muon_Pz",&T_Gen_Muon_Pz);
    chain->SetBranchAddress("T_Gen_Muon_tpPt",&T_Gen_Muon_tpPt);
    chain->SetBranchAddress("T_Gen_Muon_PDGid",&T_Gen_Muon_PDGid);
    chain->SetBranchAddress("T_Gen_Muon_status",&T_Gen_Muon_status);
    chain->SetBranchAddress("T_Gen_Muon_MotherID",&T_Gen_Muon_MotherID);
    chain->SetBranchAddress("T_Gen_Muon_FoundSTA",&T_Gen_Muon_FoundSTA);
    chain->SetBranchAddress("T_Gen_Muon_FoundL2",&T_Gen_Muon_FoundL2);
    chain->SetBranchAddress("T_Gen_Muon_L2crudeMaching",&T_Gen_Muon_L2crudeMaching);
    chain->SetBranchAddress("T_Muon_Eta",&T_Muon_Eta);
    chain->SetBranchAddress("T_Muon_Phi",&T_Muon_Phi);
    chain->SetBranchAddress("T_Muon_Pt",&T_Muon_Pt);
    chain->SetBranchAddress("T_Muon_IsGlobalMuon",&T_Muon_IsGlobalMuon);
    chain->SetBranchAddress("T_Muon_IsTrackerMuon",&T_Muon_IsTrackerMuon);
    int NbEntries = chain->GetEntries();
    cout << "nb tot events=" << NbEntries << endl;
    for (int i=0 ;i<NbEntries ; i++){
        std::vector<float> *T_Gen_Muon_tpPtCorrected  = new std::vector<float>;
        std::vector<int> *T_Gen_Muon_FoundL2Corrected = new std::vector<int>;
        std::vector<int> *T_Gen_Muon_FoundSTACorrected = new std::vector<int>;
        std::vector<int> *T_Gen_Muon_L2crudeMachingCorrected = new std::vector<int>;

        chain->GetEntry(i);
        TLorentzVector sumMuons;
        cout << "nb Gen muons =" << T_Gen_Muon_Eta->size() << endl;
        //cout << "nb of matching candidate=" << T_Gen_Muon_FoundL2->size() << endl;
        cout << "nb of tracking particles=" << T_Gen_Muon_tpPt->size() << endl;
        /// fix the L2 matching tree
        cout << "i=" << i << endl;
        int ite2 = 0;
        for (int ite = 0 ; ite < T_Gen_Muon_Eta->size() ; ite++){
            if (ite2>=T_Gen_Muon_tpPt->size()){
                T_Gen_Muon_tpPtCorrected->push_back(-1);
                T_Gen_Muon_FoundL2Corrected->push_back(0);
                T_Gen_Muon_FoundSTACorrected->push_back(0);
                T_Gen_Muon_L2crudeMachingCorrected->push_back(0);
                continue;
            }
            cout << "PTGEN=" << T_Gen_Muon_Pt->at(ite) << endl;
            float diffPt = fabs(T_Gen_Muon_Pt->at(ite)-T_Gen_Muon_tpPt->at(ite2));
            cout << "diffPt=" << diffPt << endl;
            if (diffPt>0){
                T_Gen_Muon_tpPtCorrected->push_back(-1);
                T_Gen_Muon_FoundL2Corrected->push_back(0);
                T_Gen_Muon_FoundSTACorrected->push_back(0);
                T_Gen_Muon_L2crudeMachingCorrected->push_back(0);
            }
            else{
                T_Gen_Muon_tpPtCorrected->push_back(T_Gen_Muon_tpPt->at(ite2));
                T_Gen_Muon_FoundL2Corrected->push_back(T_Gen_Muon_FoundL2->at(ite2));
                T_Gen_Muon_FoundSTACorrected->push_back(T_Gen_Muon_FoundSTA->at(ite2));
                T_Gen_Muon_L2crudeMachingCorrected->push_back(T_Gen_Muon_L2crudeMaching->at(ite2));
                ite2++;
            }
        }
   /*     for (int ite = 0 ; ite < T_Gen_Muon_Pt->size() ; ite++){
            cout << "PTtk=" << T_Gen_Muon_tpPtCorrected->at(ite) << endl;
            cout << "PTGEN=" << T_Gen_Muon_Pt->at(ite) << endl;
            cout << "L2=" << T_Gen_Muon_FoundL2Corrected->at(ite) << endl;
            cout << "STA=" << T_Gen_Muon_FoundSTACorrected->at(ite) << endl;
            cout << "L2CRUDE=" << T_Gen_Muon_L2crudeMachingCorrected->at(ite) << endl;
        }*/
        

        
        for (int j = 0  ; j < T_Gen_Muon_Eta->size() ; j++){
            if (!((fabs(T_Gen_Muon_PDGid->at(j))==13)&&(T_Gen_Muon_status->at(j)==1))) continue;
            //cout << "pt=" << T_Gen_Muon_Pt->at(j) << " eta=" << T_Gen_Muon_Eta->at(j) << " phi=" << T_Gen_Muon_Phi->at(j) << endl;
            //cout << "PDGid=" << T_Gen_Muon_PDGid->at(j) << " status=" << T_Gen_Muon_status->at(j) << " motherID=" << T_Gen_Muon_MotherID->at(j) << endl;
            TLorentzVector* muon1 = new TLorentzVector(T_Gen_Muon_Px->at(j), T_Gen_Muon_Py->at(j), T_Gen_Muon_Pz->at(j), T_Gen_Muon_Energy->at(j));
            for (int k= 0 ; k < T_Gen_Muon_Eta->size() ; k++){
                if (j==k) continue;
                if (!((fabs(T_Gen_Muon_PDGid->at(k))==13)&&(T_Gen_Muon_status->at(k)==1))) continue;
                TLorentzVector* muon2 = new TLorentzVector(T_Gen_Muon_Px->at(k), T_Gen_Muon_Py->at(k), T_Gen_Muon_Pz->at(k), T_Gen_Muon_Energy->at(k));
                sumMuons = *muon1 + *muon2;
                //cout << "Mass=" << sumMuons.M() << endl;
                //cout << "Pt=" << sumMuons.Pt() << endl;
                mu_eta = T_Gen_Muon_Eta->at(k);
                mu_pt  = T_Gen_Muon_Pt->at(k);
                pair_mass = sumMuons.M();
                pair_pt = sumMuons.Pt();
                pair_dR = deltaR(T_Gen_Muon_Phi->at(j), T_Gen_Muon_Phi->at(k), T_Gen_Muon_Eta->at(j),  T_Gen_Muon_Eta->at(k));
                //cout << "pt=" << T_Gen_Muon_Pt->at(k) << " eta=" << T_Gen_Muon_Eta->at(k) << " phi=" << T_Gen_Muon_Phi->at(k) << endl;

                foundL2 = T_Gen_Muon_FoundL2Corrected->at(k);
                foundSTA  = T_Gen_Muon_FoundSTACorrected->at(k);
                foundL2crude = T_Gen_Muon_L2crudeMachingCorrected->at(k);
                
                int nbRecoMuons = T_Muon_Pt->size();
                int iteMinDr = -1;
                float minDr = 1000;
                for (int iteReco = 0 ; iteReco < nbRecoMuons ; iteReco++){
                    float theDeltaR = deltaR(  T_Gen_Muon_Phi->at(k), T_Muon_Phi->at(iteReco), T_Gen_Muon_Eta->at(k), T_Muon_Eta->at(iteReco));
                    if (theDeltaR<minDr){
                        minDr = theDeltaR;
                        iteMinDr = iteReco;
                    }
                   // cout << " RECO pt=" << T_Muon_Pt->at(iteReco) << " eta=" << T_Muon_Eta->at(iteReco) << " phi=" << T_Muon_Phi->at(iteReco) << endl;
                }
                if (minDr < 0.1){
                    foundRECO = 1;
                    foundRECOglobal = T_Muon_IsGlobalMuon->at(iteMinDr);
                    foundRECOtracker = T_Muon_IsTrackerMuon->at(iteMinDr);
                    
                    //cout << "isGlobal=" << T_Muon_IsGlobalMuon->at(iteMinDr) << endl;
                }
                else {
                    foundRECO = 0;
                    foundRECOglobal = 0;
                    foundRECOtracker = 0;
                }
                mytree->Fill();
            }
        }
        delete T_Gen_Muon_tpPtCorrected;
        delete T_Gen_Muon_FoundL2Corrected;
        delete T_Gen_Muon_FoundSTACorrected;
        delete T_Gen_Muon_L2crudeMachingCorrected;

    }
    myOutFile->Write();
    myOutFile->Close();
    delete myOutFile;
}
