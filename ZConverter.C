//
// Created by Pranay Jain on 25/05/21.
//

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
int EM_LONG_SEG(int rodnum) {
    int seg, row = rodnum / 28;
    if (row < 8) {
        seg = 0;
    } else if (row < 17) {
        seg = 1;
    } else {
        seg = 2;
    }
    return seg;
}

// 12 radiator gaps in each hadronic module (6 per segment), and 3 hadronic modules for a total of 36 hadronic gaps (6 total segments)
// 28 rods per gap
// long_seg:row_start-row_stop = HAD1:0-5, HAD2:6-11, HAD3:12-17, HAD4:18-23, HAD5:24-29, HAD6:30-35
int HAD_LONG_SEG(int rodnum, int modnum) {
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
    vector<int> *RPD_nCherenkovs = 0;
    vector < vector < int > * > zdcRodNb;
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

    // Begin loop over events -------------------------------------------------------------------------
    for (int q = 0; q < nEntries; q++) {
        if (q % 10 == 0) cout << "\r" << left << Form("Processing event %d of %d", q, nEntries) << flush << endl;

        // Retrieve entry from trees -------------------------------------------------------------------------
        chain_event.GetEntry(q);
        chain_rpd.GetEntry(q);
        for (int k = 0; k < 4; k++) {
            ZDCchain[k]->GetEntry(q);
        }

        for (int mod = 0; mod < 4; mod++) {//start module loop
            for (int hit = 0; hit < zdcRodNb[mod]->size(); hit++) {//start hit loop

                if (mod == 0) EM_seg[EM_LONG_SEG(zdcRodNb[mod]->at(hit))]++;
                else HAD_seg[HAD_LONG_SEG(zdcRodNb[mod]->at(hit), mod)]++;

            }//end hit loop
        }//end module loop

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

// reverse commit test
