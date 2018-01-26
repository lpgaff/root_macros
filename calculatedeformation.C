// A macro to calculate the beta2, beta3 and beta4 deformation parameters

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "Riostream.h"
#endif

// Declare external variables
static double Z, A, R0;
static const int nSteps = 5000;
static int nLoops = 100;
static bool fCanvas = false;
static bool fGraphs = false;

// Define TGraphs and TCanvas
static TGraph *g2;
static TGraph *g3;
static TGraph *g4;
static TCanvas *c1;
void createCanvas() {
	c1 = new TCanvas("chisqr","chisqr",800,700);
	c1->Divide(1,3);
	fCanvas = true;
}
void createGraphs() {
	g2 = new TGraph(nSteps);
	g3 = new TGraph(nSteps);
	g4 = new TGraph(nSteps);
	g2->SetTitle("#beta_{2} vs #chi^{2};#beta_{2};#chi^{2}");
	g3->SetTitle("#beta_{3} vs #chi^{2};#beta_{3};#chi^{2}");
	g4->SetTitle("#beta_{4} vs #chi^{2};#beta_{4};#chi^{2}");
	fGraphs = true;
}

// Deformation parameters
double beta2, beta3, beta4;
double beta2bar_exp, beta3bar_exp, beta4bar_exp;

// Set parameters for grid search
float low_2, low_3, low_4;
float high_2, high_3, high_4;
float stepsize2, stepsize3, stepsize4;
int pos = 0;
int minbin;
double chisqr;
double prevchi[3] = {999.,999.,999.};


// Definition of multipole moments
double beta2bar(double beta2, double beta3, double beta4){

	double beta2bar_calc = TMath::Sqrt(5./TMath::Pi());
	beta2bar_calc *= ( (2./7.)*beta2*beta2 + (4./15.)*beta3*beta3
			+ (20./77.)*beta4*beta4 + (12./(7.*TMath::Sqrt(5.)))*beta2*beta4 );
	beta2bar_calc += beta2;
	return beta2bar_calc;

}

double beta3bar(double beta2, double beta3, double beta4){

	double beta3bar_calc = 5./TMath::Sqrt(4*TMath::Pi());
	beta3bar_calc *= ( (4./TMath::Sqrt(45.))*beta2*beta3 + (6./11.)*beta3*beta4 );
	beta3bar_calc += beta3;
	return beta3bar_calc;

}

double beta4bar(double beta2, double beta3, double beta4){

	double beta4bar_calc = 6./TMath::Sqrt(4*TMath::Pi());
	beta4bar_calc *= ( (3./7.)*beta2*beta2 + (3./11.)*beta3*beta3
			+ (243./1001.)*beta4*beta4 + (20.*TMath::Sqrt(5.)/77.)*beta2*beta4 );
	beta4bar_calc += beta4;
	return beta4bar_calc;

}

void calculatedeformation( double Zin = 88, double Ain = 224, double Q2_exp = 6.31, /*eb*/ 
				double Q3_exp = 2.57, /*eb^3/2*/ double Q4_exp = 2.60 /*eb^2*/ ){

	// Constants
	Z = Zin;
	A = Ain;
	R0 = 0.12*TMath::Power(A,1./3.);	// b^1/2
	
	// Beta bar (first order)
	beta2bar_exp = Q2_exp * TMath::Sqrt(5*TMath::Pi()) / (3*Z*TMath::Power(R0,2));
	beta3bar_exp = Q3_exp * TMath::Sqrt(7*TMath::Pi()) / (3*Z*TMath::Power(R0,3));
	beta4bar_exp = Q4_exp * TMath::Sqrt(9*TMath::Pi()) / (3*Z*TMath::Power(R0,4));
	
	// Starting values, use betaLbar
	beta2 = beta2bar_exp;
	beta3 = beta3bar_exp;
	beta4 = beta4bar_exp;
//	beta4 = 0.0803; // fix beta4 for 224Ra
	beta4 = 0.00155; // fix beta4 for 220Rn

	// Redefine step size and limits
	low_2 = beta2-0.25; high_2 = beta2+0.25;
	low_3 = beta3-0.25; high_3 = beta3+0.25;
	low_4 = beta4-0.25; high_4 = beta4+0.25;
	stepsize2 = (high_2-low_2) / nSteps;
	stepsize3 = (high_3-low_3) / nSteps;
	stepsize4 = (high_4-low_4) / nSteps;

	// Create graphs if they don't already exist
	if(!fGraphs) createGraphs();

	// Do search for best values (I think it should solve to zero)
	for(int N=1; N<=nLoops; N++){

		if(N==1){ // create canvas if we're on the first loop
			if(!fCanvas) createCanvas();
		} // if(N==1)
		
		// Vary beta2 parameter and test chisqr
		for(float i=low_2; i<=high_2; i+=stepsize2){
		
			chisqr = (beta2bar(i,beta3,beta4) - beta2bar_exp) * (beta2bar(i,beta3,beta4) - beta2bar_exp) / beta2bar_exp;
			chisqr += (beta3bar(i,beta3,beta4) - beta3bar_exp) * (beta3bar(i,beta3,beta4) - beta3bar_exp) / beta3bar_exp;
//			chisqr += (beta4bar(i,beta3,beta4) - beta4bar_exp) * (beta4bar(i,beta3,beta4) - beta4bar_exp) / beta4bar_exp;
			g2->SetPoint(pos,i,chisqr);
			pos++;
			
		} // end beta2 search
		
		// Set new beta2 to one with the best chisqr
		minbin = TMath::LocMin(g2->GetN(), g2->GetY());
		chisqr = g2->GetY()[minbin];
		if(chisqr < prevchi[0]) beta2 = g2->GetX()[minbin];
		prevchi[0] = chisqr; // reset chisqr
		pos = 0; // reset the data point number to zero
		
		// Vary beta3 parameter and test chisqr
		for(float i=low_3; i<=high_3; i+=stepsize3){
		
			chisqr = (beta2bar(beta2,i,beta4) - beta2bar_exp) * (beta2bar(beta2,i,beta4) - beta2bar_exp) / beta2bar_exp;
			chisqr += (beta3bar(beta2,i,beta4) - beta3bar_exp) * (beta3bar(beta2,i,beta4) - beta3bar_exp) / beta3bar_exp;
//			chisqr += (beta4bar(beta2,i,beta4) - beta4bar_exp) * (beta4bar(beta2,i,beta4) - beta4bar_exp) / beta4bar_exp;
			g3->SetPoint(pos,i,chisqr);
			pos++;
			
		} // end beta3 search
		
		// Set new beta3 to one with the best chisqr
		minbin = TMath::LocMin(g3->GetN(), g3->GetY());
		chisqr = g3->GetY()[minbin];
		if(chisqr < prevchi[1]) beta3 = g3->GetX()[minbin];
		prevchi[1] = chisqr; // reset chisqr
		pos = 0; // reset the data point number to zero
		
		// Vary beta4 parameter and test chisqr
		for(float i=low_4; i<=high_4; i+=stepsize4){
		
			chisqr = (beta2bar(beta2,beta3,i) - beta2bar_exp) * (beta2bar(beta2,beta3,i) - beta2bar_exp) / beta2bar_exp;
			chisqr += (beta3bar(beta2,beta3,i) - beta3bar_exp) * (beta3bar(beta2,beta3,i) - beta3bar_exp) / beta3bar_exp;
//			chisqr += (beta4bar(beta2,beta3,i) - beta4bar_exp) * (beta4bar(beta2,beta3,i) - beta4bar_exp) / beta4bar_exp;
			g4->SetPoint(pos,i,chisqr);
			pos++;
			
		} // end beta4 search
		
		// Set new beta4 to one with the best chisqr
		minbin = TMath::LocMin(g4->GetN(), g4->GetY());
//		chisqr = g4->GetY()[minbin];
//		if(chisqr < prevchi[2]) beta4 = g4->GetX()[minbin];
//		prevchi[2] = chisqr; // reset chisqr
		pos = 0; // reset the data point number to zero
				
		// Draw plots
		c1->cd(1);
		g2->Draw("ACP");
		c1->cd(2);
		g3->Draw("ACP");
		c1->cd(3);
		g4->Draw("ACP");
		
	} // for(int N=0; N<nLoops; n++)
	
	// Print out some stuff
	cout << "beta2 = " << beta2 << endl;
	cout << "beta3 = " << beta3 << endl;
	cout << "beta4 = " << beta4 << endl;
	cout << "\n              calc\t exp\t diff\n";
	cout << "beta2bar:    " << beta2bar(beta2, beta3, beta4) << "\t" << beta2bar_exp << "\t" << 100*(beta2bar(beta2, beta3, beta4)-beta2bar_exp)/beta2bar_exp << "%\n";
	cout << "beta3bar:    " << beta3bar(beta2, beta3, beta4) << "\t" << beta3bar_exp << "\t" << 100*(beta3bar(beta2, beta3, beta4)-beta3bar_exp)/beta3bar_exp << "%\n";
	cout << "beta4bar:    " << beta4bar(beta2, beta3, beta4) << "\t" << beta4bar_exp << "\t" << 100*(beta4bar(beta2, beta3, beta4)-beta4bar_exp)/beta4bar_exp << "%\n";

	return;

}