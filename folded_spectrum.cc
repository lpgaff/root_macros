/* A macro to read in relative yields and simulate a spectra */
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


double GetSigma( double FWHM ){

	return FWHM / (2.*TMath::Sqrt(2.*TMath::Log(2)));

}


void folded_spectrum( string filename, int FWHM, string outfile ){

	/* Defines */
	gStyle->SetOptLogy(1);
	double energy, yield;
	double integral = 0.;
	double w,x,y,z;
	time_t curtime;
	TRandom3 rGen(time(&curtime));
	TH1F *spec = new TH1F( "spec", ";Q-value [keV];Counts [1/keV]", 4096, -2047.5, 2048.5 );

	// Open file
	ifstream yields_in( filename.c_str(), ios::in );
	
	// Read data from file
	cout << "Energy\tYield" << endl;
	yields_in >> energy >> yield;
	while( !yields_in.eof() ) {
		
		// Loop over yield line and fill spectra
		cout << energy << "\t" << yield << endl;
		
		for( int i = 0; i < yield; i++ ){
			
			y = rGen.Gaus( energy, GetSigma( FWHM ) );
			spec->Fill(y);
			integral++;
		
		}

		yields_in >> energy >> yield;
		
	}	

	// Write to file
	string rootfilename = outfile + ".root";
	TFile *out = new TFile( rootfilename.c_str(), "RECREATE" );
	out->cd();
	spec->Write();
	yields_in.close();
	
	// Draw
	string pdffilename = outfile + ".pdf";
	gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas();
	c1->cd();
	spec->GetXaxis()->SetRangeUser(-140.,50.);
	spec->SetLineColor(kBlack);
	spec->Draw("hist");
	c1->SaveAs( pdffilename.c_str() );
	
	return;

}
