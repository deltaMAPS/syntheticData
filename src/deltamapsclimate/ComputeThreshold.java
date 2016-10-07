package deltamapsclimate;

/**
 *
 * @author Ilias.Fountalis
 */
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Set;


public class ComputeThreshold {
    
    
    double[][][] dataArray; //the 3D array (matlab output) containing grid cell time series 
    boolean maskArray[][]; //entry true if cell (i,j) is on land (or we do not want to use it)

    double alpha;
    
    int dimX,dimY,dimT;
        
    double delta;
    int nSamples = 10000;
    
    public ComputeThreshold(DeltaMapsClimate parent)
    {
        this.dataArray = parent.dataArray;
        this.maskArray = parent.maskArray;
        this.dimT = parent.dimT;
        this.dimX = parent.dimX;
        this.dimY = parent.dimY;
        this.alpha = parent.alpha;
        
        System.out.println("Setting threshold using Barlett");
        System.out.println("Using "+nSamples+" random samples.");
        this.delta = getThresholdBarlett();
                
        System.out.println("Threshold set to: "+this.delta);
    }
    
    private double getThresholdBarlett()
    {
        System.out.println("Calculating threshold using Barlett, significance level set to: "+alpha);
        int currentSamples = 0;
        double avgThreshold = 0.0, nCorrs = 0.0;
        double minSignificantCorrel = 1.0;
                    
        while(currentSamples < nSamples)
        {
            //get two grid cells
            int x1 = (int)Math.floor(Math.random()*(dimX-1));
            int x2 = (int)Math.floor(Math.random()*(dimX-1));
            int y1 = (int)Math.floor(Math.random()*(dimY-1));
            int y2 = (int)Math.floor(Math.random()*(dimY-1));
            //if these grid cells are not masked
            if(!maskArray[x1][y1] && !maskArray[x2][y2] && (x1!=x2 || y1!=y2))
            {
                if(currentSamples%1000 == 0)
                    System.out.println("Progress: "+((double)currentSamples/(double)nSamples)+"%");
                //get their cross correlation
                double[] ts1 = dataArray[x1][y1];
                double[] ts2 = dataArray[x2][y2];
                double corr = pearsonCorrel(ts1, ts2);
             
                //estimate variance using Barlett
                double var = varianceBarlett(ts1, ts2);
                //fisher transform the correlation
                //double corrF = fisherZ(corr);
                //normalize by std of Barlett
                double zijF = corr/Math.sqrt(var);
                //test significance
                double pval = 1-DistLib.normal.cumulative(zijF, 0.0, 1.0);                
                if(pval < alpha)
                {  //passed
                    avgThreshold+=corr;
                    nCorrs++;        
                    if(minSignificantCorrel > corr)
                        minSignificantCorrel = corr;
                }

                currentSamples++;
            }            
        }

        System.out.println("Min significant correlation: "+minSignificantCorrel);
        return avgThreshold/nCorrs;
    }
    
    
    private double varianceBarlett(double[] ts11, double[] ts21)
    {
        double[] ts1 = Arrays.copyOf(ts11, dimT);
        double[] ts2 = Arrays.copyOf(ts21, dimT);
        //autocorrelation from [-T,T] for ts1 and ts2
        double[] autocorrTS1 = new double[2*dimT+1];
        double[] autocorrTS2 = new double[2*dimT+1];
        for(int i = 0; i < dimT; i++)
        {
            autocorrTS1[dimT+i] = autocorr(ts1, i);
            autocorrTS2[dimT+i] = autocorr(ts2, i);
            autocorrTS1[dimT-i] = autocorrTS1[dimT+i];
            autocorrTS2[dimT-i] = autocorrTS2[dimT+i];
        }
        double var = 0.0;
        for(int i = 0; i < autocorrTS1.length; i++)
            var+= autocorrTS1[i]*autocorrTS2[i];
        return var/(double)dimT;        
    }
    
    /*
    returns the autocorrelation of ts at the specified lag
    */
    private double autocorr(double[] ts, int lag)
    {
        double mean = getMean(ts);
        double var = getVariance(ts, mean);
        double auto = 0.0;
        for(int i = 0; i < (dimT-lag); i++)        
            auto+= (ts[i]-mean)*(ts[i+lag]-mean);
        auto = auto/( ((double)dimT)*var );
        return auto;
    }
    

    private double pearsonCorrel(double[] ts1, double[] ts2)
    {
       double[] scores1 = Arrays.copyOf(ts1, dimT);
       double[] scores2 = Arrays.copyOf(ts2, dimT);
       scores1 = normalizeToZeroMeanUnitVar(scores1);
       scores2 = normalizeToZeroMeanUnitVar(scores2);
       double correl = 0.0;
       for(int i = 0; i < scores1.length; i++)
           correl+=scores1[i]*scores2[i];
       correl = correl/(double)dimT;
       
       return correl;
    }
    
    private double[] normalizeToZeroMeanUnitVar(double[] ts)
    {
        double mean = getMean(ts);
        double std = Math.sqrt(getVariance(ts, mean));
        for(int i = 0; i < ts.length; i++)
            ts[i] = (ts[i]-mean)/std;
        return ts;
    }
    
    private double getMean(double[] ts)
    {
        double mean = 0.0;
        for(int i = 0; i < ts.length; i++)
            mean+=ts[i];
        mean/=(double)(ts.length);
        return mean;
    }
       
    
    private double getVariance(double[] ts, double mean)
    {
        double var = 0.0;
        for(int i = 0; i < ts.length; i++)        
            var += Math.pow(ts[i]-mean, 2);
        var/=(double)(ts.length-1);
        return var;
    }
    
}

