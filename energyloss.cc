#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <TGaxis.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

string convertInt( int number ) {
	
	stringstream ss;
	ss << number;
	return ss.str();
	
}

string convertFloat( float number ) {
	
	stringstream ss;
	ss << number;
	return ss.str();
	
}

double GetBTh( double angle, double Ap, double At, double Eb, double Ex = 241.0 ) {

	// Returns theta angle of B using angle of T
	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1.+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);
	double x, y, BTh;
	y = TMath::Tan(angle); // y = tan(Theta_targetlab)
	x = (y*y*epsilon - TMath::Sqrt(-y*y*epsilon*epsilon + y*y + 1.) ) / (1.+y*y);
	//	cout << "Centre of mass angle: " << TMath::ACos(x)*TMath::RadToDeg() << endl;
	BTh = TMath::ATan( TMath::Sqrt( 1. - x*x ) / ( tau*epsilon + x ) );
	if( BTh < 0 ) BTh += TMath::Pi();
	//	cout << "Simulated beam angle: " << BTh*TMath::RadToDeg() << endl;
	return BTh;
	
}

double GetCOMTh( double angle, double Ap, double At, double Eb, double Ex = 241.0 ) {

	// Returns centre of mass scattering angle using angle of T
	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1.+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);
	double x, y, COMTh;
	y = TMath::Tan(angle); // y = tan(Theta_targetlab)
	x = (y*y*epsilon - TMath::Sqrt(-y*y*epsilon*epsilon + y*y + 1.) ) / (1.+y*y);
	COMTh = TMath::ACos(x);
	return COMTh;
	
}

double BTELoss_pars[6] = { -0.88, 0.91, -0.064, -27., 4.5, -0.16 };
double TTELoss_pars[6] = { -0.88, 0.91, -0.064, -27., 4.5, -0.16 };
double BSELoss_pars[6] = { 3.6, -0.15, 0.010, -37., 6.1, -0.22 };
double TSELoss_pars[6] = { 3.6, -0.15, 0.010, -37., 6.1, -0.22 };

double SP_function( double *x, double *par ) { // energy units of keV
    
    double SP_nucl = par[0] + par[1] * TMath::Log( x[0] );
    SP_nucl += par[2] * TMath::Log( x[0] ) * TMath::Log( x[0] );
    SP_nucl = TMath::Exp( SP_nucl );
    
    double SP_elec = par[3] + par[4] * TMath::Log( x[0] );
    SP_elec += par[5] * TMath::Log( x[0] ) * TMath::Log( x[0] );
    SP_elec = TMath::Exp( SP_elec );
    
    return SP_nucl + SP_elec;
    
}

bool stoppingpowers( Int_t Zp, Int_t Zt, Double_t Ap, Double_t At, string opt ) {
    
    string gElName[110] = {
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
    
    string srimfile = "../srim/"; // prefix
    if( opt.substr(0,1) == "B" ) srimfile += convertInt(Ap) + gElName[Zp-1];
    else if( opt.substr(0,1) == "T" ) srimfile += convertInt(At) + gElName[Zt-1];
    else {
        cout << "opt must equal BT, TT, BS or TS \n";
        return false;
    }
    
    if( opt.substr(1,1) == "T" ) srimfile += "_" + convertInt(At) + gElName[Zt-1] + ".txt";
    else if( opt.substr(1,1) == "S" ) srimfile += "_Si.txt";
    else {
        cout << "opt must equal BT, TT, BS or TS \n";
        return false;
    }
    
    ifstream infile;
    infile.open( srimfile.c_str(), ios::in );
    
    if( !infile.is_open() ) {
        
        cout << "Cannot open " << srimfile << endl;
        return false;
        
    }
    
    TGraph *gSP = new TGraph();
    string title = "Stopping powers for ";
    if( opt.substr(0,1) == "B" ) title += convertInt(Ap) + gElName[Zp-1];
    else if( opt.substr(0,1) == "T" )  title += convertInt(At) + gElName[Zt-1];
    if( opt.substr(1,1) == "T" ) title += " in " + convertInt(At) + gElName[Zt-1];
    else if( opt.substr(1,1) == "S" ) title += " in the Si dead layer";
    title += ";Ion energy [keV];Stopping power [MeV/(mg/cm^2)]";
    gSP->SetTitle( title.c_str() );
    
    double limits[2] = { 8E2, 1E6 };
    TF1 *fSP = new TF1( "fSP", SP_function, limits[0], limits[1], 6 );
    
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
        
        if( units == "eV" ) BEn *= 1E-3;
        else if( units == "keV" ) BEn *= 1E0;
        else if( units == "MeV" ) BEn *= 1E3;
        else if( units == "GeV" ) BEn *= 1E6;
        
        total = nucl + elec ; // MeV / ( mg / cm^2 )
        
        gSP->SetPoint( p, BEn, total );
        
        // Get next line
        getline( infile, line );
        p++;
        
        // If we've reached the end, stop
        if( line.substr( 0, 9 ) == "---------" ) endflag = true;
        if( line.substr( 0, 9 ) == " Multiply" ) endflag = true;
        
    }
    
    if( opt == "BT" ) fSP->SetParameters( BTELoss_pars );
    else if( opt == "TT" ) fSP->SetParameters( TTELoss_pars );
    else if( opt == "BS" ) fSP->SetParameters( BSELoss_pars );
    else if( opt == "TS" ) fSP->SetParameters( TSELoss_pars );
    
    gSP->Fit( fSP, "QRWM" );
    
    if( opt == "BT" ) fSP->GetParameters( BTELoss_pars );
    else if( opt == "TT" ) fSP->GetParameters( TTELoss_pars );
    else if( opt == "BS" ) fSP->GetParameters( BSELoss_pars );
    else if( opt == "TS" ) fSP->GetParameters( TSELoss_pars );
    
    
    TCanvas *c = new TCanvas();
    gSP->Draw("A*");
    gSP->GetXaxis()->SetTitleOffset(1.3);
    gSP->GetYaxis()->SetTitleOffset(1.3);
    TGaxis::SetMaxDigits(3);
    string pdfname = srimfile.substr( 0, srimfile.find_last_of(".") ) + ".pdf";
    c->SetLogx();
    c->SaveAs( pdfname.c_str() );
    
    delete c;
	
	cout << endl << convertInt(Ap) << gElName[Zp-1] << " on " << convertInt(At) << gElName[Zt-1] << endl;
    
    return true;
    
}

void energyloss( Int_t Zp = 86, Int_t Zt = 50, Double_t Ap = 220, Double_t At = 120, Double_t Eb = 2830, Double_t thick = 2.0, Double_t angle = 15.0, Int_t Nmeshpoints = 20 ) {
	
	// Fit stopping power curves
	if( !stoppingpowers( Zp, Zt, Ap, At, "BT" ) ) {
		
		cerr << "Couldn't get stopping powers\n";
		return;
		
	}

		
	// Get path length assuming interaction in centre of the target
	double dist = thick * 0.5 * ( 1. + 1./TMath::Cos( angle*TMath::DegToRad() ) );
	
	// Parameters for energy loss calculations
	double dedx = 0;
	double dx = dist/(double)(10*Nmeshpoints);
	double Ei = Eb*Ap, Ej;
	double E[1] = {Ei};
	
	// Output in gosia format
	stringstream gosia_en, gosia_dedx;
	
	// Calculate energy loss through target and get range
	for( int i = 0; i < Nmeshpoints*10; i++ ){
		
		if( E[0] < 0. ){ // if we become negative, the beam is stopped!
			
			E[0] = 0;
			break;
			
		}
		
		dedx = SP_function( E, BTELoss_pars );
		E[0] -= 1000.*dedx*dx;
		
	}
	
	// Exit energy
	Ej = E[0];
	cout << "Strip angle = " << angle << " deg\n";
	cout << "Beam angle = " << TMath::RadToDeg()*GetBTh( angle*TMath::DegToRad(), Ap, At, Eb ) << " deg\n";
	cout << "Thickness = " << thick << " mg/cm2\n";
	cout << "Effective thickness = " << dist << " mg/cm2\n";
	cout << "Average energy = " << (Ej+Ei)*0.5/1000. << " MeV\n";
	cout << "Energy range is = " << Ej/1000. << " - " << Ei/1000. << " MeV\n";

	// Calculate energy loss at each meshpoint and print results
	cout << "\nE\tdE/dx\n";
	for( int i = 0; i < Nmeshpoints; i++ ){
		
		E[0] = Ej + (Ei-Ej)*(double)i/(double)(Nmeshpoints-1);
		dedx = SP_function( E, BTELoss_pars );
		
		cout << E[0]/1000. << "\t" << dedx << endl;
		
		if( i == 0 ) gosia_en << E[0]/1000.;
		else gosia_en << ", " << E[0]/1000.;
		
		if( i == 0 ) gosia_dedx << dedx;
		else gosia_dedx << ", " << dedx;

	}

	cout << endl << gosia_en.str() << endl << gosia_dedx.str() << endl;
	
	return;
	
}