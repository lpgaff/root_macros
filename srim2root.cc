// Code to read in stopping powers from SRIM and plot in ROOT
// Liam Gaffney (liam.gaffney@cern.ch) - July 2017
//

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
using namespace std;

string convertInt( int number ) {
	
	stringstream ss;
	ss << number;
	return ss.str();
	
}

string convertFloat( float number, int precision ) {
	
	stringstream ss;
	ss << setprecision( precision ) << number;
	return ss.str();
	
}

const string gElName[110] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
    "Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",
    "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
    "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
    "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
    "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
    "Mt","Ds" };

double BELoss_pars[6] = { -5.6, 2.8, -0.2, 1.8, 0.07, -0.07 };

double SP_function( double *x, double *par ) { // energy units of keV
 
    double SP_nucl = par[0] + par[1] * TMath::Log( x[0] );
    SP_nucl += par[2] * TMath::Log( x[0] ) * TMath::Log( x[0] );
    SP_nucl = TMath::Exp( SP_nucl );
    
    double SP_elec = par[3] + par[4] * TMath::Log( x[0] );
    SP_elec += par[5] * TMath::Log( x[0] ) * TMath::Log( x[0] );
    SP_elec = TMath::Exp( SP_elec );

    return SP_nucl + SP_elec;

}

bool stoppingpowers( int Zp, int Zt, double Ap, double At, string srim_dir, TGraph *gSP, TF1 *fSP ) {
    
    string srimfile = srim_dir + "/"; // prefix
    srimfile += convertInt(Ap) + gElName[Zp-1];
    srimfile += "_" + convertInt(At) + gElName[Zt-1] + ".txt";

    ifstream infile;
    infile.open( srimfile.c_str(), ios::in );
    
    if( !infile.is_open() ) {
        
        cout << "Cannot open " << srimfile << endl;
        return false;
        
    }
    
    string title = "Stopping powers for ";
    title += convertInt(Ap) + gElName[Zp-1];
    title += " in " + convertInt(At) + gElName[Zt-1];
    title += ";Ion energy [MeV];Stopping power [MeV/(mg/cm^2)]";
    gSP->SetTitle( title.c_str() );
	
    string line, units, tmp_str;
    stringstream line_ss;
    bool endflag = false;
    double BEn, nucl, elec, total, tmp_dbl;
    int p = 0;
    
    // Test file format
    getline( infile, line );
    if( line.substr( 0, 5 ) == " ====" ) {
        
        while( line.substr( 0, 5 ) != "  ---" )
            getline( infile, line );
        
        getline( infile, line ); // read first line of data
        
    }
    
    while( !infile.eof() && !endflag ) {
        
        // Read in data
        line_ss.str("");
        line_ss << line;
        line_ss >> BEn >> units >> nucl >> elec >> tmp_dbl >> tmp_str >> tmp_dbl >> tmp_str;
        
        if( units == "eV" ) BEn *= 1E-6;
        else if( units == "keV" ) BEn *= 1E-3;
        else if( units == "MeV" ) BEn *= 1E0;
        else if( units == "GeV" ) BEn *= 1E3;
        
        total = nucl + elec ; // MeV / ( mg / cm^2 )
        
        gSP->SetPoint( p, BEn, total );
        //cout << p << "\t" << BEn << "\t" << total << endl;
        
        // Get next line
        getline( infile, line );
        p++;
        
        // If we've reached the end, stop
        if( line.substr( 0, 9 ) == "---------" ) endflag = true;
        if( line.substr( 0, 9 ) == " Multiply" ) endflag = true;
        
    }
    
    fSP->SetParameters( BELoss_pars );

    gSP->Fit( fSP, "QRWM" );

    fSP->GetParameters( BELoss_pars );
	
    TCanvas *c = new TCanvas();
    gSP->Draw("A*");
    gSP->GetXaxis()->SetTitleOffset(1.3);
    gSP->GetYaxis()->SetTitleOffset(1.3);
    TGaxis::SetMaxDigits(3);
    string pdfname = srimfile.substr( 0, srimfile.find_last_of(".") ) + ".pdf";
    c->SetLogx();
    c->SaveAs( pdfname.c_str() );
    
    return true;
    
}

void srim2root( int Zp, int Zt, double Ap, double At, string srim_dir = "../srim" ) {
		
    // Suppress some message from root
    gErrorIgnoreLevel = kWarning;
    
    // Check we have sensible elements
    if( Zp > 110 || Zt > 110 ) {
        cout << "Super heavy elements!" << endl;
        return;
    }
    
	// Open output file
	string outname = convertInt(Ap) + gElName[Zp-1] + "_" + convertInt(At) + gElName[Zt-1] + "_srim.root";
	TFile *out = new TFile(outname.c_str(),"RECREATE");

	// Graphs and functions
	TGraph *gSP = new TGraph();
	TF1 *fSP = new TF1( "fSP", SP_function, 0.5, 1000., 6 );
	
    // Setup stopping powers
    if( !stoppingpowers( Zp, Zt, Ap, At, srim_dir, gSP, fSP ) ) return;
	
	// Write plot to file
	gSP->Write("gSP");
	out->Write();

}

void energyloss( double Eb = 100., double thick = 1.0, int steps = 20 ) {
	
	
	TGraph *gSP = (TGraph*)gDirectory->Get("gSP");

	// Loop over steps through target
	double E[1] = {Eb};
	double dx = thick / (double)steps;
	double dedx;
	
	cout << "E [MeV]\tdEdx [MeV/(mg/cm2)]" << endl;
	for( int i = 0; i < steps; i++ ){
		
		dedx = gSP->Eval( E[0] );
		
		cout << E[0] << "\t" << dedx << endl;
		
		E[0] -= dedx * dx;
		
	}
	
	cout << "Initial energy = " << Eb << " MeV" << endl;
	cout << "  Final energy = " << E[0] << " MeV" << endl;
	cout << "   Energy loss = " << Eb-E[0] << " MeV" << endl;

	
}