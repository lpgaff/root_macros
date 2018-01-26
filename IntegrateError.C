// ROOT macro to integrate peaks and sum up the errors on each bin

void IntegrateError( TH1F *hist, int low_bin, int upp_bin, 
							int low_err=-1, int upp_err=-1, bool print=true ) {

	double dummyA, dummyB;
	IntegrateErrorCalc( hist, low_bin, upp_bin, dummyA, dummyB, low_err, upp_err);
	return;

}

void IntegrateError( TH1F *hist, int low_bin, int upp_bin, double &peakArea, 
							double &peakError, int low_err=-1, int upp_err=-1, bool print=true ) {

	IntegrateErrorCalc( hist, low_bin, upp_bin, peakArea, peakError, low_err, upp_err, false );
	return;

}

void IntegrateErrorCalc( TH1F *hist, int low_bin, int upp_bin, double &peakArea, 
								double &peakError, int low_err, int upp_err, bool print=true ) {

	if( upp_bin < low_bin ) {
	
		cout << "Limits wrong way around?\n";
		return;
	
	}
	
	double count = 0, error = 0;
	double bgl = 0, bgl_err = 0;
	double bgr = 0, bgr_err = 0;
	int Np = 0, Nr = 0;
	double bg_frac = 1.0;
	
	for( int i = low_bin; i <= upp_bin; i++ ) {
	
		count += hist->GetBinContent(i);
		error += hist->GetBinError(i) * hist->GetBinError(i); // sum of squares
		Np++;
	
	}
	
	if( low_err!=0 || upp_err!=0 ) { // don't count background if low_err == 0

		if( low_err < 0 ) low_err = 2*low_bin-upp_bin-1;
		if( upp_err < 0 ) upp_err = 2*upp_bin-low_bin-1;
	
		for( int j = low_err; j <= low_bin-1; j++ ) {
	
			bgl += hist->GetBinContent(j);
			bgl_err += hist->GetBinError(j) * hist->GetBinError(j); // sum of squares
			Nr++;
	
		}
	
		for( int k = upp_bin+1; k <= upp_err; k++ ) {
	
			bgr += hist->GetBinContent(k);
			bgr_err += hist->GetBinError(k) * hist->GetBinError(k); // sum of squares
			Nr++;
	
		}
	
		bg_frac = (double)Np / (double)Nr;
	
	}
	
	error = TMath::Sqrt(error);
	bgl_err = TMath::Sqrt(bgl_err);
	bgr_err = TMath::Sqrt(bgr_err);
	
	peakArea = count - bg_frac*(bgl+bgr);
	peakError = TMath::Sqrt( error*error + bg_frac*(bgl_err*bgl_err + bgr_err*bgr_err) );
	
	if( print==true ) {
	
		cout << "Central\t= " << count << " +/- " << error << endl;
		cout << "Bg L\t= " << bgl << " +/- " << bgl_err << endl;
		cout << "Bg R\t= " << bgr << " +/- " << bgr_err << endl;
		cout << "Peak\t= " << peakArea << " +/- " << peakError << endl;
	
	}
	
	return;

}