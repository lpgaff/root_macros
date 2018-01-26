// A macro to read in relative yields and simulate a spectra
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


double GetSigma( double Energy ){

	double FWHM = 2.4; // keV at 0keV
	return (0.0019*Energy+FWHM)/(2.*TMath::Sqrt(2.*TMath::Log(2)));

}

double GetEff( double *x, double *par ) {

	double E0 = 0.325;
	double N = par[5]; // normalisation
	double eff = 0;
	
	for( int i = 0; i < 5; i++ ) {

		eff += par[i] * TMath::Power( TMath::Log( x[0] / E0 ), i );
		
	}
	
	eff = N * TMath::Exp( eff );
	
	return eff;

}

void YieldSim150Nd( string filename = "150Nd_84Kr_yields.dat" ) {

	// Defines
	double energy, yield, icc;
	string mult;
	double w1,w2,w3,x,y,z;
	double prob;
	time_t curtime;
	TRandom3 rGen(time(&curtime));
	TH1F *hGammas = new TH1F( "hGammas", "84Kr on 150Nd Simulated Yields", 2000, 0, 2000 );
	TH1F *hElectrons = new TH1F( "hElectrons", "84Kr on 150Nd Simulated Yields", 2000, 0, 2000 );
	TH1F *hE0 = new TH1F( "hE0", "84Kr on 150Nd Simulated Yields", 2000, 0, 2000 );
	
	
	// Detector efficiency parameters
	double en[1];
	double ge_eff_par[6] = { 4.2383, -0.529409, -0.0321632, 0.0920407, -0.0481681, 2.32E-03 };
	double si_eff_par[6] = { 4.2383, -0.529409, -0.0321632, 0.0920407, -0.0481681, 5.00E-04 };
	TF1 *fEff = new TF1( "fEff", GetEff, 0.0, 2.0, 6 );
	fEff->SetParameters( ge_eff_par );
	TGraph *gGeEff = new TGraph( fEff );
	gGeEff->SetName( "gGeEff" );
	gGeEff->SetTitle( "Germanium efficiency;Energy [MeV];Efficiency" );
	fEff->SetParameters( si_eff_par );
	TGraph *gSiEff = new TGraph( fEff );
	gSiEff->SetName( "gSiEff" );
	gSiEff->SetTitle( "Silicon efficiency;Energy [MeV];Efficiency" );

	// Open file
	ifstream yields_in( filename.c_str(), ios::in );

	// Define x-rays
	double floP = 0.918; // flouresence probablilty
	double KbindE = 43.5689; // K binding energy
	double LbindE = 6.7215; // L binding energy (average)
	double bindE;
	
	vector< double > KxrayE, KxrayI, LxrayE, LxrayI;
	double Kprob = 0, Lprob = 0;
	KxrayE.push_back( 37.361 );		KxrayI.push_back( 47.2 );
	KxrayE.push_back( 36.847 );		KxrayI.push_back( 26.2 );
	KxrayE.push_back( 36.443 );		KxrayI.push_back( 0.00531 );
	KxrayE.push_back( 42.272 );		KxrayI.push_back( 8.97 );
	KxrayE.push_back( 43.335 );		KxrayI.push_back( 2.93 );
	KxrayE.push_back( 42.166 );		KxrayI.push_back( 4.65 );
	KxrayE.push_back( 43.451 );		KxrayI.push_back( 0.030 );
	KxrayE.push_back( 42.580 );		KxrayI.push_back( 0.121 );
	KxrayE.push_back( 43.548 );		KxrayI.push_back( 0.42 );

	LxrayE.push_back( 5.230 );		LxrayI.push_back( 5.2 );
	LxrayE.push_back( 5.208 );		LxrayI.push_back( 0.57 );
	LxrayE.push_back( 5.722 );		LxrayI.push_back( 3.12 );
	LxrayE.push_back( 6.090 );		LxrayI.push_back( 1.10 );
	LxrayE.push_back( 5.829 );		LxrayI.push_back( 0.121 );
	LxrayE.push_back( 5.723 );		LxrayI.push_back( 0.072 );
	LxrayE.push_back( 5.893 );		LxrayI.push_back( 0.0454 );
	LxrayE.push_back( 6.604 );		LxrayI.push_back( 0.50 );
	LxrayE.push_back( 6.883 );		LxrayI.push_back( 0.022 );
	LxrayE.push_back( 6.901 );		LxrayI.push_back( 0.031 );
	LxrayE.push_back( 5.146 );		LxrayI.push_back( 0.081 );
	LxrayE.push_back( 4.633 );		LxrayI.push_back( 0.213 );
	
	for( int i = 0; i < KxrayI.size(); i++ ) Kprob += KxrayI[i];
	for( int i = 0; i < LxrayI.size(); i++ ) Lprob += LxrayI[i];
	
	// Read data from file
	cout << "Energy\tYield\ticc\tmult" << endl;
	yields_in >> energy >> yield >> icc >> mult;

	while( !yields_in.eof() ) {
		
		// Loop over yield line and fill spectra
		cout << energy*1000. << "\t" << yield << "\t" << icc << "\t" << mult << endl;
		
		// gammas
		if( mult != "E0" ) {
		
			for( int i = 0; i < yield; i++ ) {
				
				w1 = rGen.Rndm(i);
				y = rGen.Gaus( energy*1000., GetSigma(energy*1000.) );
				en[0] = y/1000.;
				if( w1 < GetEff(en,ge_eff_par) ) hGammas->Fill(y);
				
			}
		
		}
		
		// electrons and x-rays
		for( int j = 0; j < yield*icc; j++ ) {
			
			w1 = rGen.Rndm() * ( Kprob + Lprob ); // K or L?

			if( mult == "E0" ) {
			
				bindE = 0;
				prob = Kprob + Lprob + 1.0;
				
			}

			else if( w1 < Kprob ) {
			
				bindE = KbindE;
				prob = 0;
				
			}
			
			else {
			
				bindE = LbindE;
				prob = Kprob;
			
			}
			
			x = energy*1000. - bindE;
			y = rGen.Gaus( x, GetSigma( x )*2.5 );
			w3 = rGen.Rndm();
			en[0] = y/1000.;
			if( w3 < GetEff(en,si_eff_par) ) {

				hElectrons->Fill( y );
				if( mult == "E0" ) hE0->Fill( y );

			}

			w2 = rGen.Rndm();
			if( w2 > floP ) continue;
			
			if( w1 < Kprob ) {
			
				for( int e = 0; e < KxrayE.size(); e++ ) {
					
					if( w1 > prob && w1 <= KxrayI[e]+prob ) {
						
						x = KxrayE[e];
						y = rGen.Gaus( x, GetSigma(x) );
						
					}
					
					prob += KxrayI[e];
					
				}
				
			}
			
			else if( w1 < Kprob + Lprob ) {
			
				for( int e = 0; e < LxrayE.size(); e++ ) {
					
					if( w1 > prob && w1 <= LxrayI[e]+prob ) {
						
						x = LxrayE[e];
						y = rGen.Gaus( x, GetSigma(x) );
						
					}
					
					prob += LxrayI[e];
					
				}
			
			}
			
			w3 = rGen.Rndm();
			en[0] = y/1000.;
			if( w3 < GetEff(en,ge_eff_par) ) hGammas->Fill(y);
			
		}

		yields_in >> energy >> yield >> icc >> mult;

	}
		
	
	// Write to file
	string outfile = filename.substr( 0, filename.find_last_of(".") ) + ".root";
	TFile *out = new TFile( outfile.c_str(), "RECREATE" );
	out->cd();
	hGammas->Write();
	hElectrons->Write();
	hE0->Write();
	gGeEff->Write();
	gSiEff->Write();
	yields_in.close();
	
	return;

}
