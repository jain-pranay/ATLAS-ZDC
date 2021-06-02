// Command Line Tool to Convert TBranches to PDF format.

#include <iostream>
#include <TTree.h>
#include <TROOT.h>
#include <vector>
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TBox.h"

string OUTPUT_FILENAME(string file_name, string branch_name) {
    file_name.erase(file_name.find_last_of("_"));
    file_name.erase(file_name.find_last_of("_"));
    return file_name + "_" + branch_name + ".pdf"; // change filename extension to change format to PDF/PNG
}

void TreeExporter() {
    // File Processing
    string filename;

    cout << "Enter name of the output file." << endl;
    cin >> filename;

    // Input TTree
    TFile *input = new TFile(Form("%s", filename.c_str()));
    TTree *tree = (TTree*)input->Get("Run4Tree");

    // Creates TCanvas to display branches.
    TCanvas *canvas = new TCanvas("c", "c", 1);

    array<string, 6> branch_names = {"EM_Row",
                                     "EM_Column",
                                     "HAD_Row",
                                     "HAD_Column",
                                     "Total_Row",
                                     "Total_Column"
                                     };

    for (auto& branch_name : branch_names) {
        tree->Draw(Form("%s", branch_name.c_str()));
        canvas->SaveAs(Form("%s", OUTPUT_FILENAME(filename, branch_name).c_str()));
    }
}