package deltamapsclimate;

import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author Ilias Fountalis.
 */
public class DeltaMapsClimate {

    String folderOut;//location of the output folder (changes dynamically every time a new data set is loaded).
    int dimX,dimY,dimT; //dimensions of .mat file (latitude/longitude/time)
    int neighSizeK;//neighborhood size for peak identification
    double alpha;//significance level for ifnering threshold delta
    double delta;// area identification homogeneity threshold
    int maxLag;// max lag for network inference
    
    int defaultMask = -1000000; //default mask value
    
    double fdrQ;//false discovery rate Q for network inference
    double[][][] dataArray;// three dimensional array where the mat file is stored
    boolean[][] maskArray;//two dimensional array (lat x lon), an entry is true if a grid cell is masked.
    
    
    public static void main(String[] args) {
        System.out.println("Input:");
        System.out.println("1. File containing a single line that has the path to the folder where the .mat files are located.");        
        System.out.println("2. Neighborhood size for peak identification (K).");
        System.out.println("3. Area identification significance level/");
        System.out.println("4. Max lag for network inference.");
        System.out.println("5. False discovery rate q");        
        System.out.println("6. Area Identification Threshold (optional)");
        
        if(args.length == 5)
        {
            String matFileLoc = args[0];
            int neighSizeK = Integer.parseInt(args[1]);
            double alpha = Double.parseDouble(args[2]);
            int maxLag = Integer.parseInt(args[3]);
            double fdrQ = Double.parseDouble(args[4]);
            DeltaMapsClimate dmc = new DeltaMapsClimate(matFileLoc, neighSizeK, alpha, maxLag, fdrQ, -1.0);
        }
        else
        {
            String matFileLoc = args[0];
            int neighSizeK = Integer.parseInt(args[1]);
            double alpha = Double.parseDouble(args[2]);
            int maxLag = Integer.parseInt(args[3]);
            double fdrQ = Double.parseDouble(args[4]);
            double delta = Double.parseDouble(args[5]);
            DeltaMapsClimate dmc = new DeltaMapsClimate(matFileLoc, neighSizeK, alpha, maxLag, fdrQ, delta);
        }                                
    }
    
    public DeltaMapsClimate(String matFileLoc, int neighSizeK, double alpha, int maxLag, double fdrQ, double delta)
    {
        long startTime = System.currentTimeMillis();        
        this.neighSizeK = neighSizeK;
        this.alpha = alpha;
        this.maxLag = maxLag;
        this.fdrQ = fdrQ;
        this.delta = delta;
        //load all the .mat file names.
        File[] matFiles = getFiles(matFileLoc);
        System.out.println("Files Loaded: "+matFiles.length);
        System.out.println();
        for(int i = 0; i < matFiles.length; i++)
        {
            //Create a new folder for the current mat file.
            File matfile = matFiles[i];            
            File dir = new File(matfile.getName());
            if(dir.mkdir())            
                this.folderOut = dir.getAbsolutePath().concat("/");            
            else 
            {
                System.out.println("Folder creation failed");
                System.exit(1);
            }
            //open the current matfile and store its contents to the dataArray
            loadFileMatlab(matfile);
            //initialize mask (default mask value is -1000000)            
            Runtime.getRuntime().gc();
            initMaskArray(defaultMask);
            //next step is to get the area identification threshold if not specified by the user.
            if(delta==-1.0)
            {
                ComputeThreshold ct = new ComputeThreshold(this);
                this.delta = ct.delta;
            }
            //peak identification
            PeakIdentification pi = new PeakIdentification(this);
            //area identification
            AreaIdentification ai = new AreaIdentification(this, pi.peaks);
            //network inference
            NetworkInferenceFDR ni = new NetworkInferenceFDR(this,ai.areas,matfile);
        }
        long endTime= System.currentTimeMillis();
        System.out.println("Finished in: "+  ((endTime-startTime)/1000)+" seconds");

    }
    
    private void loadFileMatlab(File filename)
    {        
        try
        {
            MatFileReader mat = new MatFileReader(filename);
            HashMap<String,MLArray> contents = new HashMap<String,MLArray>();
            contents = (HashMap)mat.getContent();//get all the arrays of the mat file.
            //get their names.
            Set keys = contents.keySet();
            Object[] keynames = keys.toArray();
            //Now I just want the first array.
            
            MLDouble mldouble = (MLDouble)(contents.get(keynames[0].toString()));
            double [][] temp = mldouble.getArray();
            
            int[] dimIds = mldouble.getDimensions();
            dimX = dimIds[0]; dimY = dimIds[1]; dimT = dimIds[2];
            //reshape temp
            double[][][] tempDataArray = new double[dimX][dimY][dimT];
            
            for(int i = 0; i < dimX; i++)
            {
                int timePos = 0, lonPos = 0;;
                for(int j = 0; j < temp[i].length; j++)
                {
                    tempDataArray[i][lonPos][timePos] = temp[i][j];
                    lonPos++;
                    if(lonPos%(dimY)==0) //-1 cause have counts from 0.
                    {
                        lonPos = 0;
                        timePos = timePos+1;
                    }
                }
            }
            temp = null;
            System.gc();
            dataArray = new double[dimX][dimY][dimT-2*maxLag];
            for(int i = 0; i < dimX; i++)
            {
                for(int j = 0; j < dimY; j++)
                {
                    int cc = 0;
                    for(int z = maxLag; z < dimT-maxLag; z++)
                    {
                        dataArray[i][j][cc] = tempDataArray[i][j][z];
                        cc++;
                    }
                }
            }

            dimT = dimT-2*maxLag;
            tempDataArray = null;
            System.gc();
        }
        catch(IOException io)
        {
            System.out.println("Error parsing .mat file.");
            System.out.println(io.getMessage());
            System.exit(0);
        }        
    }
    
    
    private void initMaskArray(double maskVal)
    {
        maskArray = new boolean[dimX][dimY];
        for(int i = 0; i < dimX; i++)
            for(int j = 0; j < dimY; j++)
            {
                if(dataArray[i][j][0] == maskVal)
                    maskArray[i][j] = true;
            }        
    }

    
    private File[] getFiles(String filelist)
    {
        File files[] = null;
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filelist));
            //RandomAccessFile in = new RandomAccessFile(filelist, "rw");
            File dir = new File(in.readLine());
            in.close();
            if(dir.isDirectory())
            {
                String contents[] = dir.list();
                files = new File[contents.length];
                //for(int i = 0; i < contents.length; i++)
                //    files[i] = new File(dir.getAbsolutePath()+"\\\\"+contents[i]);
                //FOR LINUX....
                for(int i = 0; i < contents.length; i++)
                    files[i] = new File(dir.getAbsolutePath()+"/"+contents[i]);
            }
        }
        catch(IOException ioe)
        {}
        return files;
    }


    
}
