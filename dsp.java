import org.jtransforms.fft.DoubleFFT_1D;

import pl.edu.icm.jlargearrays.DoubleLargeArray;


public class YShenCorrelater
{
	private static DoubleLargeArray myFFT (DoubleLargeArray input) {
		DoubleFFT_1D fftDo = new DoubleFFT_1D(input.length());
		DoubleLargeArray fft = new DoubleLargeArray(input.length() * 2);
		for(long k = 0; k < input.length(); k ++) {
			fft.setDouble(k, input.getDouble(k));
		}
        fftDo.realForwardFull(fft);
        return fft;
	}
	
	private static DoubleLargeArray myIFFT(DoubleLargeArray input) {
		DoubleFFT_1D ifftDo = new DoubleFFT_1D(input.length() / 2); 
		ifftDo.complexInverse(input, true);
		return input;
	}
	
	private static double[] complexMultiply(double r1, double i1, double r2, double i2) {
		return new double[]{r1*r2 - i1*i2, r1*i2 + r2*i1};
	}

    private static DoubleLargeArray conv(DoubleLargeArray x1, DoubleLargeArray x2) 
    {
    	long fftsize = x1.length() + x2.length() - 1;
    	DoubleLargeArray newx1 = new DoubleLargeArray(fftsize);
    	DoubleLargeArray newx2 = new DoubleLargeArray(fftsize);
    	for(long k = 0; k < fftsize; k++) {
    		double val = k >= x1.length() ?0:x1.getDouble(k);
    		newx1.setDouble(k, val);
    		val = k >= x2.length() ?0:x2.getDouble(k);
    		newx2.setDouble(k, val);
    	}
    	DoubleLargeArray fft1 = myFFT(newx1);
    	DoubleLargeArray fft2 = myFFT(newx2);
        for(long i=0; i<fftsize;i++)
        {
        	double[] tmpResult = complexMultiply(fft1.getDouble(i*2), fft1.getDouble(i*2+1), fft2.getDouble(i*2), fft2.getDouble(i*2+1));
        	fft1.setDouble(2*i, tmpResult[0]);
        	fft1.setDouble(2*i+1, tmpResult[1]);
        }
        DoubleLargeArray fullResult = myIFFT(fft1);
        DoubleLargeArray validResult = new DoubleLargeArray(fullResult.length()/2-2*x2.length()+2);
        for(long k = 0;k <validResult.length();k++){
        	validResult.set(k, fullResult.getDouble((k+x2.length()-1)*2));
        }
        return validResult;
    }
    
    private static DoubleLargeArray xcorr(DoubleLargeArray x1, DoubleLargeArray x2) {
    	double tmp = 0;
    	for(long k = 0; k < x2.length()/2; k++) {
    		tmp = x2.getDouble(k);
    		x2.setDouble(k, x2.getDouble(x2.length() - 1 - k));
    		x2.setDouble(x2.length() - 1 - k, tmp);
    	}
    	DoubleLargeArray fullCorr = conv(x1, x2);
    	return fullCorr;
    }

    public static long singalShift(double[] x1, double[] x2){
    	DoubleLargeArray a = new DoubleLargeArray(x1);
    	DoubleLargeArray b = new DoubleLargeArray(x2);
    	DoubleLargeArray c = null;
    	c = xcorr(a, b);
    	long idx = 0; 
    	double val = c.getDouble(idx);
    	for(long k =0; k<c.length();  k++) {
    		if(c.getDouble(k) > val){
    			idx = k;
    			val = c.getDouble(k);
    		}
    	}
    	return idx;
    }
}
