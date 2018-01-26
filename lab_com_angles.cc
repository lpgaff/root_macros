#include "energyloss.cc"

//void lab_com_angles( Double_t Ap = 220, Double_t At = 114, Double_t Eb = 2530.736, Double_t Ex = 241.0, Double_t low_ang = 12, Double_t upp_ang = 58, Int_t Nmeshpoints = 24 ) {
//void lab_com_angles( Double_t Ap = 220, Double_t At = 120, Double_t Eb = 2540.54, Double_t Ex = 241.0, Double_t low_ang = 12, Double_t upp_ang = 58, Int_t Nmeshpoints = 24 ) {
//void lab_com_angles( Double_t Ap = 220, Double_t At = 60, Double_t Eb = 2476.56, Double_t Ex = 241.0, Double_t low_ang = 12, Double_t upp_ang = 58, Int_t Nmeshpoints = 24 ) {
void lab_com_angles( Double_t Ap = 220, Double_t At = 60, Double_t Eb = 2476.56, Double_t Ex = 241.0, Double_t low_ang = 15, Double_t upp_ang = 52, Int_t Nmeshpoints = 37001 ) {

	
	gErrorIgnoreLevel = kWarning;
	
	double beam_ang, targ_ang, com_ang;
	double stepsize = ( upp_ang - low_ang ) / (double)( Nmeshpoints - 1 );

	cout << "Target [deg.]\tBeam [deg.]\tCoM [deg.]\n";
	for( int i = 0; i < Nmeshpoints; i++ ) {

		targ_ang = low_ang + i*stepsize;
		beam_ang = GetBTh( targ_ang*TMath::DegToRad(), Ap, At, Eb, Ex ) * TMath::RadToDeg();
		com_ang = GetCOMTh( targ_ang*TMath::DegToRad(), Ap, At, Eb, Ex ) * TMath::RadToDeg();
		
		cout << targ_ang << "\t\t" << beam_ang << "\t\t" << com_ang << endl;
	
	}
	
	return;
	
}


void cadmium_angles() {

	double beam_ang, targ_ang;
	double strip_ave;
	
	vector< double > low_strip, upp_strip;
	low_strip.push_back( 0 );		upp_strip.push_back( 1 );
	low_strip.push_back( 2 );		upp_strip.push_back( 3 );
	low_strip.push_back( 4 );		upp_strip.push_back( 5 );
	low_strip.push_back( 6 );		upp_strip.push_back( 7 );
	low_strip.push_back( 8 );		upp_strip.push_back( 9 );
	low_strip.push_back( 10 );		upp_strip.push_back( 12 );

	cout << "Exp\tTarget [deg.]\tBeam ave.\tTarg. low\tTarg. upp\n";
	for( int i = 0; i < low_strip.size(); i++ ) {

		strip_ave = 0.5*( low_strip[i] + upp_strip[i] );
		
		// Average
		targ_ang = TMath::RadToDeg() * TMath::ATan((9.+(15.5-strip_ave)*2.)/32.2);
		beam_ang = GetBTh( targ_ang*TMath::DegToRad(), 220., 114., 563.611, 241. ) * TMath::RadToDeg();
		cout << i << "\t" << targ_ang << "\t\t" << beam_ang << "\t\t";

		targ_ang = TMath::RadToDeg() * TMath::ATan((9.+(15.05-upp_strip[i])*2.)/32.2);
		cout << targ_ang << "\t\t";

		targ_ang = TMath::RadToDeg() * TMath::ATan((9.+(15.95-low_strip[i])*2.)/32.2);
		cout << targ_ang << "\n";

	}

	return;

}