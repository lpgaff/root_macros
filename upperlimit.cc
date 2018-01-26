// Liam Gaffney (Liam.Gaffney@fys.kuleuven.be) 17/11/2014
// Quick ROOT script to solve the upper limit according to O. Helene, NIM 212 (1983) 319-322
// Equations and method stolen from Simon Sels' powerpoint presentation (16/05/2013)


#include "TMath.h"
#include "TF1.h"

#define GAUSLIM 20.0 // limit (infinity!) for the gaussian-like function, sigma = 1
#define NSTEPS 1000 // number of steps to search for A

// alpha = I[(A-a)/sigma] / I[-a/sigma]
// or
// I[(A-a)/sigma] = alpha * I[-a/sigma]
// where
// I(z) = ( 1 / 2*pi ) * int^{inf}_{z} exp(-0.5*x^2) dx
 
void upperlimit( double a, double sigma, double alpha = 0.317310508 ) {

	TF1 *g = new TF1( "g", "TMath::Exp(-0.5*x*x)/TMath::Sqrt(2.0*TMath::Pi())", -GAUSLIM, GAUSLIM );
	
	double A;	// upper limit, to search
	double maxA = 2.*a + 10.*sigma;	// guess the maximum
	double stepA = maxA / (double)NSTEPS;
	
	double LHS;
	double RHS = alpha * g->Integral( -a / sigma, GAUSLIM );
	
	cout << endl << "alpha * I[-a/sigma] = " << RHS << endl << endl;
	
	// Rough search
	cout << "Rough search; step size in A = " << stepA << endl;
	for ( int i = 1; i <= NSTEPS; i++ ) {
	
		A = i * stepA;
		LHS = g->Integral( ( A - a ) / sigma, GAUSLIM );
		
		if ( LHS < RHS ) break;
	
	}
	
	cout << "\n\tcurrent A = " << A << endl << endl;
	cout << "I[(A-a)/sigma] = " << LHS << endl << endl;
	
	// Fine search
	stepA *= (4./(double)NSTEPS); // redefine step size
	double minA = A - 0.5*(double)NSTEPS*stepA; // redefine range

	cout << "Fine search; step size in A = " << stepA << endl;	
	for ( int i = 1; i <= NSTEPS; i++ ) {
	
		A = i * stepA + minA;
		LHS = g->Integral( ( A - a ) / sigma, GAUSLIM );
		
		if ( LHS < RHS ) break;
	
	}
	
	cout << "\n\tfinal A = " << A << endl << endl;
	cout << "I[(A-a)/sigma] = " << LHS << endl << endl;
	
	return;
	
}