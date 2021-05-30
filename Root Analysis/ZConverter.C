//
// Created by Pranay Jain on 25/05/21.
//
// EM Segmentation Test

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
    int rowNum, RODS_PER_ROW = 10;
    rowNum = rodNum / RODS_PER_ROW;
    return rowNum;
}

void ZConverter() {

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
    vector<int> *RPD_nCherenkovs = 0;
    vector<vector<int>*> zdcRodNb;
    vector<int> rodNum;
    zdcRodNb.resize(4);

    int trackID;
    int EM_rows[26];
    int HAD_rows[10];
    double EM_seg[3];
    double HAD_seg[6];
    double energy;


    // Output Tree Branches
    // Vector Branches
    tOut->Branch("lastStepZ", &LastStepInVolume);
    tOut->Branch("rpdNcherenkov", &RPD_nCherenkovs);

    // Standard Branches
    tOut->Branch("energy", &energy, "energy/D");
    tOut->Branch("trackID", &trackID, "trackID/I");

    // Array Branches
    tOut->Branch("EM_seg", EM_seg, "EM_seg[3]/D");
    tOut->Branch("HAD_seg", HAD_seg, "HAD_seg[6]/D");
    tOut->Branch("EM_rows", EM_rows, "EM_rows[26]/I");
    tOut->Branch("HAD_rows", HAD_rows, "HAD_rows[10]/I");

    // Setup reading of input tree --------------------------------

    //EventData contains information related to the primary particle
    TChain chain_event("EventData");
    //RPD1tree contains nCherenkovs, a 256 length vector corresponding to the # of photons in each rod for an event
    TChain chain_rpd("RPD1tree");
    //there are 4 zdc trees corresponding to the EM (ZDC1) + 3HAD (ZDC2-4) modules
    std::vector < TChain * > ZDCchain;
    for (int i = 0; i < 4; i++) {
        ZDCchain.push_back(new TChain(Form("ZDC%dtree", i + 1)));
    }

    if (filename.find("_") != string::npos) filename.erase(filename.find_last_of("_"));

    //Add all the files to our TChains
    for (int i = 0; i < num_files; i++) {

        chain_event.Add(Form("%s_%d.root", filename.c_str(), i));
        chain_rpd.Add(Form("%s_%d.root", filename.c_str(), i));

        for (int k = 0; k < 4; k++) {
            ZDCchain[k]->Add(Form("%s_%d.root", filename.c_str(), i));
        }
    }

    //Make sure arrays are zeroed
    for (int i = 0; i < 6; i++) {
        if (i < 3) EM_seg[i] = 0;
        HAD_seg[i] = 0;
    }

    //Set addresses for tree variables
    chain_event.SetBranchAddress("lastStepZ", &LastStepInVolume);
//		chain_event.SetBranchAddress(	"energy",&energy);
    chain_rpd.SetBranchAddress("nCherenkovs", &RPD_nCherenkovs);

    for (int k = 0; k < 4; k++) {
        ZDCchain[k]->SetBranchAddress("rodNo", &zdcRodNb[k]);
    }

    int nEntries = ZDCchain[0]->GetEntries();

    // Begin loop over events
    for (int q = 0; q < nEntries; q++) {
        if (q % 10 == 0) cout << "\r" << left << Form("Processing event %d of %d", q, nEntries) << flush << endl;

        // Retrieve entry from trees
        chain_event.GetEntry(q);
        chain_rpd.GetEntry(q);
        for (int k = 0; k < 4; k++) {
            ZDCchain[k]->GetEntry(q);
        }

        for (int mod = 0; mod < 4; mod++) {//start module loop
            for (int hit = 0; hit < zdcRodNb[mod]->size(); hit++) {//start hit loop
                if (mod == 0) {
                    // EM Z segmentation
                    int hit_row = EM_Z_SEG(zdcRodNb[mod]->at(hit));
                    EM_rows[hit_row]++;
                } else {
                    // HAD Z segmentation
                }

            }//end hit loop
        }//end module loop

        int sumRows = 0;
        for (int i = -1; i < 27; i++) {
            cout << "Row " << i << ": " << EM_rows[i] << endl;
            sumRows += EM_rows[i];
        }
        cout << sumRows << endl;

        trackID = q;

        tOut->Fill();

        //zero arrays
        for (int i = 0; i < 6; i++) {
            if (i < 3) EM_seg[i] = 0;
            HAD_seg[i] = 0;
        }

    }//end event loop
    fOut->Write();
}