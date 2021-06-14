// Tree Converter File for Test Beam 2021

#include <iostream>
#include <TTree.h>
#include <TROOT.h>
#include <vector>
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TBox.h"

using namespace std;

// Segments any given rodNum into row number within the EM Module
// RODS_PER_ROW = number of rods in a single row. For RUN4 config, 29 rods.
int EM_Z_SEG(int rodNum) {
    int rowNum, RODS_PER_ROW = 29;
    rowNum = rodNum / RODS_PER_ROW;
    return rowNum;
}

// Segments any given rodNum into row number within the HAD Modules
// RODS_PER_ROW = number of rods in a single row. For RUN4 config, 29 rods.
int HAD_Z_SEG(int rodNum) {
    int rowNum, RODS_PER_ROW = 59;
    rowNum = rodNum / RODS_PER_ROW;
    return rowNum;
}

// Segments any given rodNum into column number with the EM and HAD Modules (works with both)
// RODS_PER_ROW = number of rods in a single row. For RUN4 config, 29 rods.
int EM_X_SEG(int rodNum) {
    int columnNum, RODS_PER_ROW = 29;
    columnNum = rodNum % RODS_PER_ROW;
    return columnNum;
}

int HAD_X_SEG(int rodNum) {
    int columnNum, RODS_PER_ROW = 59;
    columnNum = rodNum % RODS_PER_ROW;
    return columnNum;
}

void TestBeamTreeConverter() {

    // File Processing
    // Enter full filename, including .root, followed by a space and then number of consecutive files.
    string filename;
    int num_files;

    cout << "Enter Filename, Number of Consecutive Files" << endl;
    cin >> filename >> num_files;


    if (filename.find(".") != string::npos) filename.erase(filename.find_last_of("."));
    if (filename.find(".") != string::npos) filename.erase(filename.find_last_of("."));

    TFile *fOut = new TFile(Form("%s_TB21Out.root", filename.c_str()), "RECREATE");
    TTree *tOut = new TTree("TestBeamTree", "TestBeamTree");

    // Output Tree Variables
    vector<double> *LastStepInVolume = 0;
    vector<int>	*RPD_nCherenkovs = 0;
    vector<int> *EM_nCherenkovs = 0;
    vector<int> *HAD_nCherenkovs = 0;
    vector<vector<int>*> zdcRodNb(4);

    TH1I *EM_Row = new TH1I("EM_Row", "EM_Row", 11, 0, 11);
    TH1I *EM_Column = new TH1I("EM_Column", "EM_Column", 29, 0, 29);
    TH1I *HAD_Row = new TH1I("HAD_Row", "HAD_Row", 36, 0, 36);
    TH1I *HAD_Column = new TH1I("HAD_Column", "HAD_Column", 59, 0, 59);
    TH1I *Total_Row = new TH1I("Total_Row", "Total_Row", 47, 0, 47);
    TH1I *Total_Column = new TH1I("Total_Column", "Total_Column", 59, 0, 59);

    TH2I *EM_Cone = new TH2I("EM_Cone", "EM_Cone", 11, 0, 11, 29, 0, 29);
    TH2I *HAD_Cone = new TH2I("HAD_Cone", "HAD_Cone", 47, 0, 47, 59, 0, 59);
    gStyle->SetPalette(kRainBow);

    int trackID = 0;

    // Output Tree Branches
    // Vector Branches
    tOut->Branch("LastStepZ", &LastStepInVolume);
    tOut->Branch("RPD_nCherenkovs", &RPD_nCherenkovs);
    tOut->Branch("EM_nCherenkovs", &EM_nCherenkovs);
    tOut->Branch("HAD_nCherenkovs", &HAD_nCherenkovs);

    tOut->Branch("EM_Row", &EM_Row);
    tOut->Branch("HAD_Row", &HAD_Row);
    tOut->Branch("EM_Column", &EM_Column);
    tOut->Branch("HAD_Column", &HAD_Column);
    tOut->Branch("Total_Row", &Total_Row);
    tOut->Branch("Total_Column", &Total_Column);
    tOut->Branch("EM_Cone", &EM_Cone);
    tOut->Branch("HAD_Cone", &HAD_Cone);

    // Standard Branches
    tOut->Branch("TrackID", &trackID, "TrackID/I");


    // Setting up and reading  input tree
    // EventData contains information related to the primary particle
    TChain Event_Chain("EventData");
    // RPD1tree contains nCherenkovs, a 256 length vector corresponding to the # of photons in each rod for an event
    TChain RPD_Chain("RPD1tree");
    // There are 4 zdc trees corresponding to the EM (ZDC1) + 3HAD (ZDC2-4) modules
    vector<TChain*> ZDC_Chain;
    for (int i = 0; i < 4; i++) {
        ZDC_Chain.push_back(new TChain(Form("ZDC%dtree", i+1)));
    }

    if (filename.find("_") != string::npos) {
        filename.erase(filename.find_last_of("_") );
    }

    // Add all files to TChains
    for (int i = 0; i < num_files; i++) {
        Event_Chain.Add(Form("%s_%d.root", filename.c_str(), i)); // Adds all EventData branches
        RPD_Chain.Add(Form("%s_%d.root", filename.c_str(), i)); // Adds all RPD1tree branches
        // Loops over ZDC1tree, ZDC2tree and ZDC3tree, adds all branches
        for(int k = 0; k < 4; k++) {
            ZDC_Chain[k]->Add(Form("%s_%d.root", filename.c_str() ,i));
        }
    }

    //Set addresses for tree variables
    Event_Chain.SetBranchAddress("lastStepZ", &LastStepInVolume);
    RPD_Chain.SetBranchAddress("nCherenkovs", &RPD_nCherenkovs);
    for (int i = 0; i < 4; i++) {
        ZDC_Chain[i]->SetBranchAddress("rodNo", &zdcRodNb[i]);
        if (i == 0) {
            ZDC_Chain[i]->SetBranchAddress("nCherenkovs", &EM_nCherenkovs);
        } else {
            ZDC_Chain[i]->SetBranchAddress("nCherenkovs", &HAD_nCherenkovs);
        }
    }

    int nEntries = ZDC_Chain[0]->GetEntries();

    // Loop over number of entries (events set per job)
    for (int q = 0; q < nEntries; q++) {
        cout << "\r" << left << Form("Processing event %d of %d", q, nEntries) << flush << endl;

        // Retrieve entry from trees
        Event_Chain.GetEntry(q);
        RPD_Chain.GetEntry(q);
        for (int i = 0; i < 4; i++){
            ZDC_Chain[i]->GetEntry(q);
        }

        // Module Loop for EM + HAD1,2,3 modules
        for (int mod = 0; mod < 4; mod++){
            for (int hit = 0; hit < zdcRodNb[mod]->size(); hit++){
                // EM Module Processing
                if (mod == 0) {
                    EM_Row->Fill(EM_Z_SEG(zdcRodNb[mod]->at(hit)));
                    EM_Column->Fill(EM_X_SEG(zdcRodNb[mod]->at(hit)));
                    Total_Row->Fill(EM_Z_SEG(zdcRodNb[mod]->at(hit)));
                    Total_Column->Fill(EM_X_SEG(zdcRodNb[mod]->at(hit)) + 15);

                    EM_Cone->Fill(EM_Z_SEG(zdcRodNb[mod]->at(hit)), EM_X_SEG(zdcRodNb[mod]->at(hit)));
                }
                    // HAD Modules Processing
                else {
                    if (mod == 1) {
                        HAD_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit))); // Indexes row 0-11 as HAD1 module.
                        HAD_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 11);
                        Total_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));

                        HAD_Cone->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)), HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                    if (mod == 2) {
                        HAD_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 12); // Indexes row 12-23 as HAD2 module.
                        HAD_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 12 + 11);
                        Total_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));

                        HAD_Cone->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 12, HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                    if (mod == 3) {
                        HAD_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 24); // Indexes row 24-35 as HAD3 module.
                        HAD_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 24 + 11);
                        Total_Column->Fill(HAD_X_SEG(zdcRodNb[mod]->at(hit)));

                        HAD_Cone->Fill(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 24, HAD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                }
            }
        }
        trackID = q;
        tOut->Fill();
    }
    fOut->Write();
}
