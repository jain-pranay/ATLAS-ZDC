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

// 25 radiator gaps in an EM module (cannot divide evenly by 3 longitudinal segments)
// 28 rods per gap
// long_seg:row_start-row_stop = EM1:0-7, EM2:8-16, EM3:17-24 (arbitrary selection)
int EM_LONG_SEG(int rodnum){
    int seg, row = rodnum / 29;
    if (row < 8) {
        seg = 0;
    } else if (row < 17) {
        seg = 1;
    } else {
        seg = 2;
    }
    return seg;
}

// 12 radiator gaps in each hadron module (6 per segment), and 3 hadron modules for a total of 36 hadron gaps (6 total segments)
// 28 rods per gap
// long_seg:row_start-row_stop = HAD1:0-5, HAD2:6-11, HAD3:12-17, HAD4:18-23, HAD5:24-29, HAD6:30-35
int HAD_LONG_SEG(int rodnum, int modnum){
    int seg, row = rodnum / 29;
    seg = row / 6;
    return seg + ((modnum - 1) * 2);
}

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
    int rowNum, RODS_PER_ROW = 29;
    rowNum = rodNum / RODS_PER_ROW;
    return rowNum;
}

// Segments any given rodNum into column number with the EM and HAD Modules (works with both)
// RODS_PER_ROW = number of rods in a single row. For RUN4 config, 29 rods.
int MOD_X_SEG(int rodNum) {
    int columnNum, RODS_PER_ROW = 29;
    columnNum = rodNum % RODS_PER_ROW;
    return columnNum;
}

void Run4TreeConverter() {

    // File Processing
    // Enter full filename, including .root, followed by a space and then number of consecutive files.
    string filename;
    int num_files;

    cout << "Enter Filename, Number of Consecutive Files" << endl;
    cin >> filename >> num_files;


    if (filename.find(".") != string::npos) filename.erase(filename.find_last_of("."));
    if (filename.find(".") != string::npos) filename.erase(filename.find_last_of("."));

    TFile *fOut = new TFile(Form("%s_Out.root", filename.c_str()), "RECREATE");
    TTree *tOut = new TTree("Run4Tree", "Run4Tree");

    // Output Tree Variables
    vector<double> *LastStepInVolume = 0;
    vector<int>	*RPD_nCherenkovs = 0;
    vector<int> *EM_nCherenkovs = 0;
    vector<int> *HAD_nCherenkovs = 0;
    vector<vector<int>*> zdcRodNb;
    vector<int> *EM_Row = 0;
    vector<int> *HAD_Row = 0;
    vector<int> *EM_Column = 0;
    vector<int> *HAD_Column = 0;
    vector<int> *Total_Row = 0;
    vector<int> *Total_Column = 0;

    vector<int> rodNum;
    zdcRodNb.resize(4);

    int trackID = 0;
    double EM_Seg[3] = {0};
    double HAD_Seg[6] = {0};
    double energy = 0;


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

    // Standard Branches
    tOut->Branch("Energy", &energy, "Energy/D");
    tOut->Branch("TrackID", &trackID, "TrackID/I");

    // Array Branches
    tOut->Branch("EM_Seg", EM_Seg, "EM_Seg[3]/D");
    tOut->Branch("HAD_Seg", HAD_Seg, "HAD_Seg[6]/D");

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

    // Ensure all arrays are zeroed
    for (int i = 0; i < 6; i++) {
        if (i < 3) {
            EM_Seg[i] = 0;
        }
        HAD_Seg[i] = 0;
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
                    EM_Seg[EM_LONG_SEG(zdcRodNb[mod]->at(hit))]++;
                    EM_Row->push_back(EM_Z_SEG(zdcRodNb[mod]->at(hit)));
                    EM_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                    Total_Row->push_back(EM_Z_SEG(zdcRodNb[mod]->at(hit)));
                    Total_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                }
                // HAD Modules Processing
                else {
                    HAD_Seg[HAD_LONG_SEG(zdcRodNb[mod]->at(hit), mod)]++;
                    if (mod == 1) {
                        HAD_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit))); // Indexes row 0-11 as HAD1 module.
                        HAD_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 26);
                        Total_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                    if (mod == 2) {
                        HAD_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 12); // Indexes row 12-23 as HAD2 module.
                        HAD_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 12 + 26);
                        Total_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                    if (mod == 3) {
                        HAD_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 24); // Indexes row 24-35 as HAD3 module.
                        HAD_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                        Total_Row->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)) + 24 + 26);
                        Total_Column->push_back(MOD_X_SEG(zdcRodNb[mod]->at(hit)));
                    }
                }
            }
        }
        trackID = q;
        tOut->Fill();

        // Zero data structures again for next iteration
        for (int i = 0; i < 6; i++) {
            if (i < 3) {
                EM_Seg[i] = 0;
            }
            HAD_Seg[i] = 0;
        }
        EM_Row->clear();
        HAD_Row->clear();
        EM_Column->clear();
        HAD_Column->clear();
        Total_Row->clear();
        Total_Column->clear();
    }
    fOut->Write();
}