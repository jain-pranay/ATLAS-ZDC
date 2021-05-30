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
    int seg, row = rodnum/28;
    if (row < 8) {
        seg = 0;
    } else if (row<17) {
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
    int seg, row = rodnum / 28;
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
    int rowNum, RODS_PER_ROW = 10;
    rowNum = rodNum / RODS_PER_ROW;
    return rowNum;
}

void Run4TreeConverter() {

    // File Processing
    // Enter full filename, including .root, followed by a space and then number of consecutive files.
    string input;

    cout << "Enter Filename, Number of Consecutive Files" << endl;
    getline(cin, input);
    string filename = input.substr(0, input.find_first_of(" "));
    int num_files = stoi(input.substr(input.find_first_of(" ") + 1));

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
    vector<int> *rowem = 0;
    vector<int> *hadem = 0;
    vector<int> rodNum;
    zdcRodNb.resize(4);

    int trackID = 0;
    int EM_Rows[26] = {0};
    int HAD_Rows[10] = {0};
    double EM_Seg[3] = {0};
    double HAD_Seg[6] = {0};
    double energy = 0;


    // Output Tree Branches
    // Vector Branches
    tOut->Branch("LastStepZ", &LastStepInVolume);
    tOut->Branch("RPD_nCherenkovs", &RPD_nCherenkovs);
    tOut->Branch("EM_nCherenkovs", &EM_nCherenkovs);
    tOut->Branch("HAD_nCherenkovs", &HAD_nCherenkovs);
    tOut->Branch("rowem", &rowem);
    tOut->Branch("hadem", &hadem);

    // Standard Branches
    tOut->Branch("Energy", &energy, "Energy/D");
    tOut->Branch("TrackID", &trackID, "TrackID/I");

    // Array Branches
    tOut->Branch("EM_Seg", EM_Seg, "EM_Seg[3]/D");
    tOut->Branch("HAD_Seg", HAD_Seg, "HAD_Seg[6]/D");
    tOut->Branch("EM_Rows", EM_Rows, "EM_Rows[26]/I");
    tOut->Branch("HAD_Rows", HAD_Rows, "HAD_Rows[10]/I");

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
    for (int i = 0; i < 26; i++) {
        if (i < 10) {
            if (i < 6) {
                if (i < 3) {
                    EM_Seg[i] = 0;
                }
                HAD_Seg[i] = 0;
            }
            HAD_Rows[i] = 0;
        }
        EM_Rows[i] = 0;
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
                if (mod == 0) {
                    EM_Seg[EM_LONG_SEG(zdcRodNb[mod]->at(hit))]++;
                    EM_Rows[EM_Z_SEG(zdcRodNb[mod]->at(hit))]++;
                    rowem->push_back(EM_Z_SEG(zdcRodNb[mod]->at(hit)));
                }
                else {
                    HAD_Seg[HAD_LONG_SEG(zdcRodNb[mod]->at(hit), mod)]++;
                    HAD_Rows[HAD_Z_SEG(zdcRodNb[mod]->at(hit))]++;
                    hadem->push_back(HAD_Z_SEG(zdcRodNb[mod]->at(hit)));
                }
            }
        }
        trackID = q;
        tOut->Fill();

        // Zero arrays again for next loop
        for (int i = 0; i < 26; i++) {
            if (i < 10) {
                if (i < 6) {
                    if (i < 3) {
                        EM_Seg[i] = 0;
                    }
                    HAD_Seg[i] = 0;
                }
                HAD_Rows[i] = 0;
            }
            EM_Rows[i] = 0;
        }
        for (int i =0; i < rowem->size(); i++) {
            rowem->at(i) = 0;
        }

    }
    fOut->Write();
}