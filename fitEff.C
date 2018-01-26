// A ROOT macro to fit a Miniball efficiency curve and calculate the correct correlated uncertainty function
// Liam Gaffney April 2013 (with edits in September 2014)
// Liam.Gaffney@fys.kuleuven.be

#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#define Npar 5	// number of parameters in efficiency function = polynomial degree + 1
#define E0 325.0	// denominator in efficiency function

string convertInt( int number ) {

	stringstream ss;
	ss << number;
	return ss.str();

}

string convertDouble( double number ) {

	stringstream ss;
	ss << number;
	return ss.str();

}

// Efficiency curve function, exponential of polynomial of logs to degree Npar
Double_t fEffCurve( Double_t *x, Double_t *p ) {

	Double_t f = p[0];
	for( int i = 1; i < Npar; i++ )
		f += p[i] * TMath::Log(x[0]/E0)**i;
	f = TMath::Exp(f);
	
	return f;

}

// Total error function, including correlations
Double_t fEffCurveErr( Double_t *x, Double_t *p ) {

	double effpar[Npar];
	for( int m = 0; m < Npar; m++ )
		effpar[m] = p[Npar*Npar+m];

	Double_t f = 0, g = 0;
	Double_t h = fEffCurve( x, effpar )**2;
	for( int m = 0; m < Npar; m++ ) {
		for( int n = 0; n < Npar; n++ ) {
			g = h*TMath::Log(x[0]/E0)**m;
			g *= TMath::Log(x[0]/E0)**n;
			g *= p[m*Npar+n];
			f += g;
		}
	}
	
	return TMath::Sqrt(f);

}

// Total error function, without correlations
Double_t fEffCurveErrNoCov( Double_t *x, Double_t *p ) {

	double effpar[Npar];
	for( int m = 0; m < Npar; m++ )
		effpar[m] = p[Npar+m];

	Double_t f = 0, g = 0;
	Double_t h = fEffCurve( x, effpar );
	for( int m = 0; m < Npar; m++ ) {
		g = h*TMath::Log(x[0]/E0)**m;
		g *= p[m];
		f += g**2;
	}
	
	return TMath::Sqrt(f);

}

void fitEff( string inputfile = "geEff_Jun2009.dat", int format = 0 ) {

	// format = 1; raw intensities
	// format = 2; absolute efficiencies
	if ( format !=1 && format !=2 ) {
	
		cout << "Incorrect format specifier:\n";
		cout << "format = 1; raw intensities\n";
		cout << "format = 2; absolute efficiencies\n";
		return;
		
	}
		
	// Open input file and store
	ifstream in;
	in.open(inputfile.c_str(),ios::in);
	if( !in.is_open() ) {
	
		cout << "Unable to open " << inputfile << endl;
		return;

	}

	// Open root output files 
	string rootname = inputfile.substr(0, inputfile.find_last_of(".")) + ".root";
	TFile *root = new TFile(rootname.c_str(), "RECREATE");
	
	// Variables
	string line;
	stringstream data (stringstream::in | stringstream::out);
	double En, dEn, Igam, dIgam, Iexp, dIexp; 
	double Eff, dEff;
	int i=0;

	// Read experimental data and fill TGraph
	TGraphErrors *gEff = new TGraphErrors();
	gEff->SetTitle("Miniball efficiency curve June 2009;Energy [keV];Efficiency [arb. units]");

	while( getline(in, line) ){
	
		if( line.substr(0,1) == '#' ) continue;
		
		data.clear();
		data << line;
		
		if ( format == 1 ) {
			data >> En >> dEn >> Igam >> dIgam >> Iexp >> dIexp;
			Eff = Iexp*0.00005 / Igam;
			dEff = Eff * TMath::Sqrt( (dIexp/Iexp)**2 + (dIgam/Igam)**2 );
		}
		else if ( format == 2 ) {
			data >> En >> dEn >> Eff >> dEff;
		}
		
		gEff->SetPoint( i, En, Eff );
		gEff->SetPointError( i, dEn, dEff );
		
		i++;
		
	}
	in.close();
	
	// Efficiency curve and partial (partial) derivative definitions
	int lowEn = 20; // upper limit on energy. Used for ranges etc
	int uppEn = 2000; // upper limit on energy. Used for ranges etc
	TF1 *effCurve = new TF1("fEffCurve", fEffCurve, (double)lowEn, (double)uppEn, Npar); // in keV
	TF1 *effCurveErr = new TF1("fEffCurveErr", fEffCurveErr, (double)lowEn, (double)uppEn, Npar*(Npar+1));
	TF1 *effCurveErrNoCov = new TF1("fEffCurveErrNoCov", fEffCurveErrNoCov, (double)lowEn, (double)uppEn, Npar*2);
	
	// Set initial parameters
	effCurve->SetParameters( 2.5, -0.7, -0.05, 0.16, -0.07 );
	effCurve->SetParNames( "a", "b", "c", "d", "e" );
	
	// Perform the fit
	TFitResultPtr res = gEff->Fit( "fEffCurve", "SRE"); // Use E for better errors
	int fitstatus = res;
	if( fitstatus < 0 ) {
	
		cout << "The fit is not OK... Try changing the initial parameters\n";
		return;
		
	}
	
	// Normalise errors to chi2 / NDF = 1
	res->NormalizeErrors();

	// Extract statistics
	double chisq = res->Chi2() / res->Ndf();
	TMatrixDSym cov = res->GetCovarianceMatrix();
	
	// Feedback parameters to error functions
	const int Nelem = cov.GetNoElements();
	double parArray[Nelem+Npar];
	double parArrayNoCov[Npar*2];
	cov.GetMatrix2Array(parArray);	
	for( int i = 0; i < Npar; i++ ) {
		parArray[Nelem+i] = res->Value(i);
		parArrayNoCov[i] = res->ParError(i);		
		parArrayNoCov[Npar+i] = res->Value(i);		
	}
	effCurveErr->SetParameters( parArray );
	effCurveErrNoCov->SetParameters( parArrayNoCov );
	
	// Save the errors and lower and upper limits to graphs
	TGraph *gEffErr = new TGraph(uppEn-lowEn);
	TGraph *gEffErrLow = new TGraph(uppEn-lowEn);
	TGraph *gEffErrUpp = new TGraph(uppEn-lowEn);
	TGraph *gEffErrNoCov = new TGraph(uppEn-lowEn);
	gEffErr->SetTitle("Miniball efficiency curve June 2009;Energy [keV];Fractional Error");
	gEffErrLow->SetTitle("Miniball efficiency curve June 2009;Energy [keV];Lower efficiency limit [arb. units]");
	gEffErrUpp->SetTitle("Miniball efficiency curve June 2009;Energy [keV];Upper efficiency limit [arb. units]");
	gEffErrNoCov->SetTitle("Miniball efficiency curve June 2009;Energy [keV];Fractional Error without covariance terms");

	for( int i = lowEn; i < uppEn; i++ ) {
		double dum, err;
		double eff = effCurve->Eval(i);
		double err = effCurveErr->Eval(i);
		gEffErr->SetPoint( i-lowEn, i, err/eff );
		gEffErrLow->SetPoint( i-lowEn, i, eff-err );
		gEffErrUpp->SetPoint( i-lowEn, i, eff+err );
		gEffErrNoCov->SetPoint( i-lowEn, i, effCurveErrNoCov->Eval(i)/eff );
	}
	
	// Compute the confidence intervals at the x points in a new graph
	TGraphErrors *gEffConf = new TGraphErrors(uppEn-lowEn);
	gEffConf->SetTitle("Miniball efficiency curve June 2009: Confidence intervals;Energy [keV];Efficiency [arb. units]");
	for( int j=lowEn; j<uppEn; j++ ) gEffConf->SetPoint(j-lowEn, j, effCurve->Eval(j));
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gEffConf);
	
	// Draw
	TCanvas *c1 = new TCanvas();
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Miniball efficiency curve June 2009");
	mg->Add(gEff,"P");
	mg->Add(gEffErrLow,"C");	// error function with covariance
	mg->Add(gEffErrUpp,"C");	// upper and lower limits
//	mg->Add(gEffConf,"4");		// confidence intervals
	mg->Draw("A");				// experimental data + errors and confidence intervals
	effCurve->Draw("CSAME");	// fitted curve
	
	// Line colours and styles
	gStyle->SetOptFit(1111);
	gEffErrLow->SetLineStyle(2);
	gEffErrUpp->SetLineStyle(2);
	gEffConf->SetFillStyle(3010);
	gEffConf->SetFillColor(kBlue);
	
	// Print and write to file
	res->Print("V");
	cout << "\nChisq / NDF = " << chisq << endl;
	string pdfname = inputfile.substr(0, inputfile.find_last_of(".")) + ".pdf";
	c1->SaveAs(pdfname.c_str());
	gEff->Write("gEff");
	gEffErr->Write("gEffErr");
	gEffErrLow->Write("gEffErrLow");
	gEffErrUpp->Write("gEffErrUpp");
	gEffErrNoCov->Write("gEffErrNoCov");
	effCurve->Write("effCurve");
	effCurveErr->Write("effCurveErr");
	res->Write();
	
	// Influence of covariance terms on the error...
	mg->Delete();
	c1->Clear();
	TMultiGraph *errcompare = new TMultiGraph();
	errcompare->Add(gEffErr);
	errcompare->Add(gEffErrNoCov);
	errcompare->SetTitle(gEffErr->GetTitle());
	errcompare->Draw("ac");
	gEffErr->SetLineColor(kRed);
	gEffErrNoCov->SetLineColor(kBlue);
	
	return;

}
