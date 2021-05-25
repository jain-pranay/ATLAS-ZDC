#include "TFile.h"
#include "TTreeReader.h"
#include <TTree.h>
#include <TROOT.h>
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TBox.h"
#include <vector>

using namespace std;

//25 radiator gaps in an EM module (cannot divide evenly by 3 longitudinal segments)
//28 rods per gap
//long_seg:row_start-row_stop = EM1:0-7, EM2:8-16, EM3:17-24 (arbitrary selection)
int EM_LONG_SEG(int rodnum){
    int seg=0;
    int row = rodnum/28;
    if(row<8)				seg = 0;
    else if(row<17)	seg = 1;
    else 						seg = 2;
    return seg;
}

//12 radiator gaps in each hadronic module (6 per segment), and 3 hadronic modules for a total of 36 hadronic gaps (6 total segments)
//28 rods per gap
//long_seg:row_start-row_stop = HAD1:0-5, HAD2:6-11, HAD3:12-17, HAD4:18-23, HAD5:24-29, HAD6:30-35
int HAD_LONG_SEG(int rodnum, int modnum){
    int seg=0;
    int row = rodnum/28;
    seg = row/6;
    return seg+((modnum-1)*2);
}

void TreeConverter() {

    //NAME OF FIRST INPUT FILE
    string filename = "run4_1p_380G_0.root";

    //NUMBER OF CONSECUTIVE FILES BEING CONVERTED
    int num_files = 1;

    //OUTPUT FILE/TREE CREATION
    string outputName = filename;
    if(outputName.find("/") != string::npos) outputName.erase( 0, outputName.find_last_of("/") + 1 );
    if(outputName.find(".") != string::npos) outputName.erase( outputName.find_last_of(".") );

    TFile *fOut = new TFile( Form("%s_Out.root",outputName.c_str()), "RECREATE" );
    TTree *tOut = new TTree("TestBeam_Tree","TestBeam_Tree");

    //TREE VARIABLES
    std::vector<double> 	*LastStepInVolume=0;
    std::vector<int>			*RPD_nCherenkovs=0;
    std::vector< std::vector<int> *>	zdcRodNb;
    zdcRodNb.resize(4);

    int 		trackID;
    double EM_seg[3];
    double HAD_seg[6];
    double energy;


    //OUTPUT BRANCH CREATION
    //VECTOR BRANCHES
    tOut->Branch("lastStepZ",  				&LastStepInVolume);
    tOut->Branch("rpdNcherenkov",  		&RPD_nCherenkovs);
    //STANDARD BRANCHES
    tOut->Branch("energy",   					&energy, 								"energy/D");
    tOut->Branch("trackID",    				&trackID,   						"trackID/I");
    //ARRAY BRANCHES
    tOut->Branch("EM_seg", 						EM_seg, 								"EM_seg[3]/D");
    tOut->Branch("HAD_seg", 					HAD_seg, 								"HAD_seg[6]/D");

    // Setup reading of input tree --------------------------------

    //EventData contains information related to the primary particle
    TChain chain_event("EventData");
    //RPD1tree contains nCherenkovs, a 256 length vector corresponding to the # of photons in each rod for an event
    TChain chain_rpd("RPD1tree");
    //there are 4 zdc trees corresponding to the EM (ZDC1) + 3HAD (ZDC2-4) modules
    std::vector<TChain*> ZDCchain;
    for(int i=0; i<4; i++ ){
        ZDCchain.push_back( new TChain( Form ("ZDC%dtree", i+1) ) ) ;
    }

    if(outputName.find("_") != string::npos) outputName.erase( outputName.find_last_of("_") );

    //Add all the files to our TChains
    for(int i = 0; i < num_files; i++){

        chain_event.Add( Form("%s_%d.root",outputName.c_str(),i) );
        chain_rpd.Add( Form("%s_%d.root",outputName.c_str(),i) );

        for(int k=0; k<4; k++ ){
            ZDCchain[k]->Add( Form("%s_%d.root",outputName.c_str(),i) );
        }
    }

    //Make sure arrays are zeroed
    for (int i=0; i < 6 ; i++){
        if(i<3) EM_seg[i] 	  = 0;
        HAD_seg[i] 						= 0 ;
    }

    //Set addresses for tree variables
    chain_event.SetBranchAddress(	"lastStepZ",&LastStepInVolume);
//		chain_event.SetBranchAddress(	"energy",&energy);
    chain_rpd.SetBranchAddress(		"nCherenkovs",&RPD_nCherenkovs);

    for(int k=0; k<4; k++ ){
        ZDCchain[k]->SetBranchAddress("rodNo", &zdcRodNb[k]);
    }

    int nEntries = ZDCchain[0]->GetEntries();

    // Begin loop over events -------------------------------------------------------------------------
    for (int q=0;q<nEntries; q++)
    {
        if(q%10 == 0) cout << "\r" << left << Form("Processing event %d of %d", q, nEntries) << flush;

        // Retrieve entry from trees -------------------------------------------------------------------------
        chain_event.GetEntry(q);
        chain_rpd.GetEntry(q);
        for(int k = 0; k < 4; k++){
            ZDCchain[k]->GetEntry(q);
        }

        for( int mod = 0; mod < 4; mod++){//start module loop
            for (int hit=0; hit < zdcRodNb[mod]->size() ; hit++){//start hit loop

                if(mod==0) 	EM_seg[EM_LONG_SEG(zdcRodNb[mod]->at(hit))]++;
                else 				HAD_seg[HAD_LONG_SEG(zdcRodNb[mod]->at(hit), mod)]++;

            }//end hit loop
        }//end module loop

        trackID		= q;

        tOut->Fill();

        //zero arrays
        for (int i=0; i < 6 ; i++){
            if(i<3) EM_seg[i] 	  = 0;
            HAD_seg[i] 						= 0;
        }

    }//end event loop

    fOut->Write();

}
