// Get B(E2) ratios in the SU(3) limit: B(E2; J+2 -> J) / B(E2; 2 -> 0)

#include "TMath.h"

double clebsch(int J /* 1 < J < 19 */){
	double coeff[19] = { 1.0,			// Ji = 0; Jf = 2;
						 0.774597,		// Ji = 1; Jf = 3;
						 0.717137,		// Ji = 2; Jf = 4;
						 0.690066,		// Ji = 3; Jf = 5;
						 0.674200,		// Ji = 4; Jf = 6;
						 0.663747,		// Ji = 5; Jf = 7;
						 0.656330,		// Ji = 6; Jf = 8;
						 0.650791,		// Ji = 7; Jf = 9;
						 0.646496,		// Ji = 8; Jf = 10;
						 0.643066,		// Ji = 9; Jf = 11;
						 0.640264,		// Ji = 10; Jf = 12;
						 0.637931,		// Ji = 11; Jf = 13;
						 0.635959,		// Ji = 12; Jf = 14;
						 0.634270,		// Ji = 13; Jf = 15;
						 0.632807,		// Ji = 14; Jf = 16;
						 0.631527,		// Ji = 15; Jf = 17;
						 0.630399,		// Ji = 16; Jf = 18;
						 0.629396,		// Ji = 17; Jf = 19;
						 0.628499};		// Ji = 18; Jf = 20;
	return coeff[J];
}

double be2ratios_su3(double J, double N, double norm){ 
	return norm * ( (J+2)*(J+1) / ( (2*J+3)*(2*J+5) ) ) * ( (2*N-J)*(2*N+J+3)*5 / ( N*(2*N+3) ) );
}

double be2ratios_clebsch(double J, double Q){ 
	/* Returns B(E2; J+2 -> J) */
	return ( (2*J+1) / (2*J+5) ) * Q*Q * clebsch(J)*clebsch(J) * 5 / (16*TMath::Pi());
}

