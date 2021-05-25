#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TGraph.h"
#include <vector>
#include <TTree.h>
#include <TROOT.h>

//HISTOGRAM VARIABLES
string axis_title[10];
string hist_title[10];

//use params to alter titles and other TH1D options
TH1D* make1DHist( double h_low, double h_high, int h_bins, int param, int unique_id ){

    TH1D* h1 = new TH1D( Form("hist_name_%d_%d", param, unique_id) ,
                         Form("%s %d", hist_title[param].c_str(), unique_id) ,
                         h_bins, h_low, h_high);
    h1->GetXaxis()->SetTitle(axis_title[param].c_str());
    h1->GetYaxis()->SetTitle(axis_title[param+1].c_str());
    h1->SetTitleSize(0.045, "XY");
    h1->GetXaxis()->SetTitleOffset(1.05);
    h1->GetYaxis()->SetTitleOffset(1.05);
    h1->GetXaxis()->SetTitleSize(0.07);
    h1->GetYaxis()->SetTitleSize(0.07);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.06);
    return h1;
}

//Refer to TH1D class for options regarding hist creation and drawing
void Draw1DPlot(vector < TH1D* > h1, TCanvas * canv1){
    canv1->cd(0);

    int markerstyle[4];
    markerstyle[0]=22;markerstyle[1]=23;markerstyle[2]=33;markerstyle[3]=8;
    int markercolor[4];
    markercolor[0]=38;markercolor[1]=6;markercolor[2]=30;markercolor[3]=46;
    string title_plot[4];
    title_plot[0]="EM Longitudinal 1";title_plot[1]="EM Longitudinal 2";
    title_plot[2]="EM Longitudinal 3";title_plot[3]="EM Total";

    auto legend10 = new TLegend(0.7,0.68,0.90,0.90);

    for(int i=0;i<4;i++){
        h1[i]->SetMarkerStyle(markerstyle[i]);
        h1[i]->SetMarkerSize(2);
        h1[i]->SetLineColorAlpha(markercolor[i],0.5);
        h1[i]->SetLineWidth(3);

        h1[i]->GetXaxis()->SetTitleSize(0.07);
        h1[i]->GetXaxis()->SetLabelSize(0.055);
        h1[i]->GetYaxis()->SetLabelSize(0.065);

        h1.at(i)->Draw("HIST SAME");
        h1.at(i)->GetYaxis()->SetRangeUser(1, 1e4);

        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.25);

        gStyle->SetTitleFontSize(0.07);
        gStyle->SetOptStat(0);

        gPad->SetLogy();
        gPad->Update();

        legend10->AddEntry(h1.at(i),Form(title_plot[i].c_str(),i),"l");
    }
    legend10->Draw();
}

void SimpleTreeReader() {
    //NAME OF INPUT FILE
    string filename = "./run4_3kg_100Gev.root";
    //READ INPUT FILE
    TFile *myFile = new TFile(filename.c_str(),"READ");
    //GET TREE FROM INPUT FILE
    TTree *TestBeam_Tree = (TTree*)myFile->Get("TestBeam_Tree");
    //NAME OF OUTPUT FILE
    TFile *fOut = new TFile( "output_file.root" , "RECREATE" );

    // Setup reading of trees ----------------------------------

    //LONGITUDINAL MODULE SEGMENTS
    int num_EM_seg = 3;
    int num_HAD_seg = 6;

    //STORE TREE DATA
    double 	energy	= 0, rpdNcheren = 0;
    std::vector<double> *ZFirstInt=0;
    std::vector<int> 		*RPD_nCherenkovs=0;
    //THREE LONGITUDINAL SEGMENTS IN EM
    double  EM_seg[num_EM_seg];
    //2 LONGITUDINAL SEGMENTS IN EACH OF 3 HADs
    double	HAD_seg[num_HAD_seg];

    //EVENT LOOP VARIABLES
    double 	total_EM = 0, total_HAD = 0, total_RPD = 0, total_light = 0;
    //HISTOGRAMS
    vector < TH1D* > h_EMlight, h_HADlight;
    vector < TH2D* > h_EMlight2D, h_HADlight2D;

    //ASSIGN BRANCHES TO CORRESPONDING VARIABLE
    TestBeam_Tree->SetBranchAddress("energy",&energy);
    TestBeam_Tree->SetBranchAddress("lastStepZ",&ZFirstInt);
    TestBeam_Tree->SetBranchAddress("rpdNcherenkov",&RPD_nCherenkovs);
    TestBeam_Tree->SetBranchAddress("EM_seg",&EM_seg);
    TestBeam_Tree->SetBranchAddress("HAD_seg",&HAD_seg);

    //PARAMS FOR HIST MAKING
    axis_title[0]="axis1"; axis_title[1]="axis2"; axis_title[2]="axis3";
    hist_title[0]="title1"; hist_title[1]="title2"; hist_title[2]="title3";

    //CREATE TH1Ds USING make1DHist FUNCTION
    for(int i=0;i<num_HAD_seg+1;i++)	{
        if(i<num_EM_seg+1) h_EMlight.push_back( make1DHist(0,1,20,0,i) );//low, high, bins, param, idx
        h_HADlight.push_back( make1DHist(0,0.05,40,1,i) );//low, high, bins, param, idx
    }

    int nEntries = TestBeam_Tree->GetEntries();

    // Begin loop over events -------------------------------------------------------------------------
    for (int q=0; q<nEntries ; q++){
        if(q%100 == 0) cout << "\r" << left << Form("Processing event %d of %d", q, nEntries) << flush;

        //RETRIEVE DATA FROM TREE
        TestBeam_Tree->GetEntry(q);

        //LOOP OVER RPD FIBERS (total_RPD not currently used, shown to inform)
        for (int fib_idx=0; fib_idx < RPD_nCherenkovs->size(); fib_idx++){
            total_RPD+=RPD_nCherenkovs->at(fib_idx);
        }

        //LOOP OVER EM/HAD SEGMENTS TO CALCULATE TOTALS
        for(int i=0;i<num_HAD_seg;i++)	{
            if(i<num_EM_seg) total_EM += EM_seg[i];
            total_HAD	+= HAD_seg[i];
        }

        //FILL HISTOGRAMS
        for(int i=0;i<num_HAD_seg;i++)	{
            if(i<num_EM_seg) h_EMlight[i]->Fill(EM_seg[i]/(total_RPD+total_EM+total_HAD));
            h_HADlight[i]->Fill(HAD_seg[i]/(total_RPD+total_EM+total_HAD));
        }

        //FILL TOTAL HISTOGRAMS
        h_EMlight[3]->Fill((total_EM)/(total_EM+total_HAD));
        h_HADlight[6]->Fill((total_HAD)/(total_EM+total_HAD));

        //Energy of Primary can be accessed using "energy" variable
        //std::cout << Form("Primary Energy #%d = %.2f [GeV] " , q, energy/1000) << std::endl;

        //ZERO EVENT VARIABLES
        total_EM 	= 0;
        total_RPD = 0;
        total_HAD = 0;

    }//End event loop
    //-------------------------------------------------------------------------

    TestBeam_Tree	->ResetBranchAddresses();

    //CREATE CANVAS AND DRAW HISTOGRAMS
    TCanvas *c_stackEMlight       = new TCanvas("c_stackEMlight","c_stackEMlight",1200,1100);
    Draw1DPlot(h_EMlight, c_stackEMlight);

    fOut->Write();
}
