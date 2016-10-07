package deltamapsclimate;

import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class NetworkInference {
    
    Vector<Area> areas;
    double[][][] dataArray;
    boolean[][] maskArray;
    double significanceLevel;    
    int dimX,dimY,dimT;
    int maxLag;
    
    double latitudes90[] = {-89,-87,-85,-83,-81,-79,-77,-75,-73,-71,-69,-67,-65,-63,-61,-59,-57,-55,-53,-51,-49,-47,-45,-43,-41,-39,-37,-35,-33,-31,-29,-27,-25,-23,-21,-19,-17,-15,-13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89};    
    double[] currentLatitudes;
    
    String folderOut;
    NewPublicFunctions pf = new NewPublicFunctions();
    
    public NetworkInference(DeltaMapsClimate parent, Vector<Area> areas, File matFile)
    {        
        this.areas = areas;
        this.significanceLevel = parent.fdrQ;
        this.maxLag = parent.maxLag;
        this.folderOut = parent.folderOut;
        this.maskArray = parent.maskArray;
        System.out.println("--------- Network Inference Starts --------------");
        
        //Step 1. load the mat file again.
        System.out.println("Loading .mat file again.");
        loadFileMatlab(matFile);
        
        System.out.println("Weighting time series in respect to grid cell size");
        setCurrentLatitudes();
        resizeTimeSeries();
        setTrueAreaSize();
        
        System.out.println("Constructing area cummulative anomaly time series.");
        constructAreaCA();
        
        System.out.println("Constructing network");
        networkInference();
        
        System.out.println("Calculating area strengths");
        setAreaStrength();                        
        
        setAreaCentroids();        
        System.out.println("Exporting maps.");
        //exportLagkMap(areas.firstElement(),"full");
        
        //COMMENT 1. 
        /*
        if you need to get all the areas then change (copy paste) the for loop to
        for(int x = 0; x < areas.size(); x++)
        */
        for(int x = 0; x < 5; x++)
        {
            exportLinkMap(areas.get(x), Integer.toString(areas.get(x).areaID));
            exportLagInfo(areas.get(x), Integer.toString(areas.get(x).areaID));
        }
        
        exportGraphCLM(areas, "");
        exportAreaIDAndCentroid();
        exportGraphNet(areas,"orig");
        
        mapColoring(areas);
        printAreaMap();
        exportForARI();
        
        //finally sort the areas by size.
        ComparatorAreaTrueSizeBarlett mycomp = new ComparatorAreaTrueSizeBarlett();
        java.util.Collections.sort(areas,mycomp);
        //get the largest area (which is typically ENSO).
        Area enso = areas.firstElement();        
        exportLinkMap(areas.firstElement(), "ENSO");
        exportLagInfo(areas.firstElement(), "ENSO");
        //print a map of the area (to make sure that it is ENSO)
        printArea(enso);
        //print the strength and size of the area
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("ensoStats.txt",true));//set true to append
            out.println(enso.trueSize+"\t"+enso.size+"\t"+enso.strength);
            out.close();
        }
        catch(IOException ioe)
        {}

        printStrengthMap();
    }
    
    private void setAreaStrength()
    {
        HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
        for(int i = 0; i < areas.size(); i++)
            edges.addAll(areas.get(i).edges);
        for(int i = 0; i < areas.size(); i++)
        {
            Area area = areas.get(i);
            double strength = 0.0;
            java.util.Iterator<AreaEdge> it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                if(e.undirected)
                {
                    if(e.from == area.areaID)
                        strength+=Math.abs(e.weight);
                }
                else
                {
                    if(e.from == area.areaID)
                        strength+=Math.abs(e.weight);
                    else if(e.to == area.areaID)
                        strength+=Math.abs(e.weight);
                }
            }
            areas.get(i).strength = strength;
        }
    }

    
    private void networkInference()
    {
        System.out.println("Network Inference Started.");
        System.out.println("Significance level: "+this.significanceLevel);
        System.out.println("Max lag: "+this.maxLag);
        int nAreas = areas.size();
        int nEdgesTotal = 0, nEdgesDirected = 0;
        for(int i = 0; i < nAreas; i++) //search for an edge between all pairs of areas.
        {
            Area alpha = areas.get(i);
            for(int j = (i+1); j < nAreas; j++)
            {
                Area beta = areas.get(j);            
                
                //step one is to get the correlogram
                double[] correlogram = getCorrelogram(alpha.areaCA, beta.areaCA);
                //this array will keep the variance estimated by Bartlett.
                double[] bvar = new double[correlogram.length];
                //step two is to identify significant correlations, use the barlett formula
                double[] significantCorrelations = new double[correlogram.length];
                //(1) get the central part of the time series
                double ts1[] = new double[alpha.areaCA.length-2*maxLag];
                double ts2[] = new double[beta.areaCA.length-2*maxLag];
    
                for(int x = maxLag; x < maxLag+ts1.length; x++)
                {
                    ts1[x-maxLag] = alpha.areaCA[x];
                    ts2[x-maxLag] = beta.areaCA[x];
                }
                //estimate the variance according to Bartlett
                int T = ts1.length;
                double[] autocorrTS1 = new double[2*T+1];
                double[] autocorrTS2 = new double[2*T+1];
                for(int x = 0; x < T; x++)
                {
                    autocorrTS1[T+x] = autocorr(ts1, x);
                    autocorrTS2[T+x] = autocorr(ts2, x);
                    autocorrTS1[T-x] = autocorrTS1[T+x];
                    autocorrTS2[T-x] = autocorrTS2[T+x];
                }
                double var = 0; //after this point var is NOT normalized by T - \tau. 
                for(int x = 0; x < autocorrTS1.length; x++)
                    var+=autocorrTS1[x]*autocorrTS2[x];
                //and now find the significant correlations
                for(int x = 0; x < correlogram.length; x++)
                {
                    int currentLag = Math.abs(maxLag-x);
                    double corr = Math.abs(correlogram[x]); //by taking the absolute corr value is like doing a two-tailed t-test                    
                    //normalize with Bartlett variance
                    corr = corr/Math.sqrt(var/(T-currentLag)); //note: here we further normalize by T-\tau.
                    bvar[x] = Math.sqrt(var/(T-currentLag));
                    double pval = 1-DistLib.normal.cumulative(corr, 0.0, 1.0);                
                    if(pval < significanceLevel)
                        significantCorrelations[x] = correlogram[x];//passed the test                                        
                }
                //finally get the edge.
                AreaEdge e = edgeInference(significantCorrelations);
                AreaEdge foo = edgeInfernce(significantCorrelations, bvar);
               
                if(e!=null)
                {   
                   
                    e.correlogram = new double[correlogram.length];
                    e.correlogram = correlogram;
                    e.bvar = bvar;
                    e.maxLag = foo.maxLag;
                    e.minLag = foo.minLag;
                    nEdgesTotal++;
                    //then we have an edge, we need to get the direction and set the weight.
                    double alphaSTD = Math.sqrt(getVariance(ts1, getMean(ts1)));
                    double betaSTD = Math.sqrt(getVariance(ts2, getMean(ts2)));
                    e.weight = alphaSTD*betaSTD*e.sij;
   
                    if(e.lag > 0 ) ///lags are positive alpha -> beta
                    {
                        e.from = alpha.areaID;
                        e.to = beta.areaID;
                        areas.get(i).edges.add(e);
                        nEdgesDirected++;
                    }
                    else if(e.lag < 0) //negative lag beta -> alpha
                    {
                        e.from = beta.areaID;
                        e.to = alpha.areaID;
                        areas.get(j).edges.add(e);
                        nEdgesDirected++;
                    }
                    else//undirected (lag range includes 0)
                    {
                        AreaEdge e2 = new AreaEdge();
                        e2.correlogram = new double[correlogram.length];
                        e2.correlogram = correlogram;
                        e2.sigcorrs = new double[correlogram.length];
                        e2.sigcorrs = e.sigcorrs;
                        e2.bvar = e.bvar;
                        e2.sij = e.sij;
                        e2.lag = e.lag;
                        e2.maxLag = e.maxLag;
                        e2.minLag = e.minLag;
                        e2.weight = e.weight;
                        e.from = alpha.areaID;
                        e.to = beta.areaID;
                        e.undirected=true;
                        areas.get(i).edges.add(e);
                        e2.from = beta.areaID;
                        e2.to = alpha.areaID;
                        e2.undirected=true;
                        areas.get(j).edges.add(e2);
                    }
                }
            }
        }
        System.out.println("#Edges (total): "+nEdgesTotal);
        System.out.println("#Edges (directed): "+nEdgesDirected);             
    }
    
    private AreaEdge edgeInfernce(double[] significantCorrelations, double[] bartlettVar)
    {
        double maxCorr  = 0;
        int lag = 0;
        for(int i = 0; i < significantCorrelations.length; i++)
            if(Math.abs(maxCorr) < Math.abs(significantCorrelations[i]))
                maxCorr = significantCorrelations[i];
        if(maxCorr == 0)
        {
            return null;
        }
        else
        {
            AreaEdge e = new AreaEdge();
            e.sij = maxCorr;
            e.lag = lag-this.maxLag;
            e.sigcorrs = new double[significantCorrelations.length];
            e.sigcorrs = significantCorrelations;
            //now get the Bartlett variance at maxCorr
            double maxBvar = bartlettVar[lag];
            //scan through the significant correlations
            int[] lagRange = new int[bartlettVar.length];
            for(int i = 0; i < significantCorrelations.length; i++)
            {
                double corr = Math.abs(significantCorrelations[i]);
                if(corr > 0)
                {
                    double cVar = bartlettVar[i];
                    if(corr >= Math.abs(maxCorr))
                    {
                        if(Math.abs(maxCorr)+maxBvar > corr-cVar)
                            lagRange[i] = 1;
                    }
                    else
                    {
                    if(Math.abs(maxCorr)-maxBvar < corr+cVar)
                        lagRange[i] = 1;
                    }
                }
            }
            int minLag = 0, maxLag = 0;
            for(int i = 0; i < lagRange.length; i++)
                if(lagRange[i] > 0)
                {
                    minLag = i;
                    break;
                }
            for(int i = lagRange.length-1; i >= 0; i--)
                if(lagRange[i] > 0)
                {
                    maxLag = i;
                    break;
                }
            e.maxLag = maxLag-this.maxLag;
            e.minLag = minLag-this.maxLag;
            return e;
        }
    }
        
    
    private AreaEdge edgeInference(double[] significantCorrelations)
    {
        double maxCorr = 0;
        int lag = 0;
        //step 1. find the maximum significant correlation in absolute sense
        for(int i = 0; i < significantCorrelations.length; i++)
            if(Math.abs(significantCorrelations[i]) > Math.abs(maxCorr))
            {
                maxCorr = significantCorrelations[i];
                lag = i;
            }

        if(maxCorr == 0)//no edge
        {
            return null;
        }
        else
        {
            AreaEdge e = new AreaEdge();
            e.sij = maxCorr;
            e.lag = lag-this.maxLag;
            e.sigcorrs = new double[significantCorrelations.length];
            e.sigcorrs = significantCorrelations;
         
            return e;
        }        
    }
    
         /*
    returns the autocorrelation of ts at the specified lag
    */
    private double autocorr(double[] tsFoo, int lag)
    {
        double[] ts = Arrays.copyOf(tsFoo, tsFoo.length);
        double mean = getMean(ts);
        double var = getVariance(ts, mean);
        double auto = 0.0;
        for(int i = 0; i < (ts.length-lag); i++)        
            auto+= (ts[i]-mean)*(ts[i+lag]-mean);
        auto = auto/( ((double)ts.length)*var );
        return auto;
    }

    private double[] getCorrelogram(double[] ts1Foo, double[] ts2Foo)
    {                        
        double[] ts1 = Arrays.copyOf(ts1Foo, ts1Foo.length);
        double[] ts2 = Arrays.copyOf(ts2Foo, ts2Foo.length);
        double cijtauArray[] = new double[2*this.maxLag+1];        
     
        //corelation at lag zero
        cijtauArray[this.maxLag]  = pearsonCorrel(ts1, ts2, maxLag,0);        
                
        //start by shifting the second time series to the right
        for(int i = 1; i <= this.maxLag; i++)                
            cijtauArray[this.maxLag+i] = pearsonCorrel(ts1, ts2,maxLag,i);            
        //continue by shifting the time serises to the left
        for(int i = -1; i >= -this.maxLag; i--)                    
            cijtauArray[this.maxLag+i] =  pearsonCorrel(ts1, ts2,maxLag,i);                    
        return cijtauArray;
    }
    
    private double pearsonCorrel(double[] scores1, double[] scores2, int startPos, int lag)
    {    
        double correl = 0.0;
    
        double[] ts1 = new double[scores1.length-2*startPos];
        double[] ts2 = new double[scores1.length-2*startPos];
        
        //get the part of the time series that you want
        for(int i = startPos; i < startPos+ts1.length; i++)
        {
            ts1[i-startPos] = scores1[i];
            ts2[i-startPos] = scores2[i+lag];
        }
        
        //set to unit variance, zero mean first
        ts1 = normalizeToZeroMeanUnitVar(ts1);
        ts2 = normalizeToZeroMeanUnitVar(ts2);
        for(int i = 0; i < ts1.length; i++)
            correl+=ts1[i]*ts2[i];
        
        return correl/(double)ts1.length;        
    }     

    
    private double[] normalizeToZeroMeanUnitVar(double[] ts)
    {
        double mean = getMean(ts);
        double std = Math.sqrt(getVariance(ts, mean));
        for(int i = 0; i < ts.length; i++)
            ts[i] = (ts[i]-mean)/std;
        return ts;
    }
    
    private double getVariance(double[] ts, double mean)
    {
        double var = 0.0;
        for(int i = 0; i < ts.length; i++)        
            var += Math.pow(ts[i]-mean, 2);
        var/=(double)(ts.length-1);
        return var;
    }

    private double getMean(double[] ts)
    {
        double mean = 0.0;
        for(int i = 0; i < ts.length; i++)
            mean+=ts[i];
        mean/=(double)(ts.length);
        return mean;
    }

    
    private void constructAreaCA()
    {
        for(int i = 0; i < areas.size(); i++)            
        {
            Area alpha = areas.get(i);
            double areaCA[] = new double[dimT];
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int cellID = alpha.cells.get(j);
                int pos[] = pf.getPosition(cellID, dimX, dimY);
                double[] ts = dataArray[pos[0]][pos[1]];
                
                for(int x = 0; x < ts.length; x++)
                    areaCA[x] = areaCA[x]+ts[x];
                
            }
            areas.get(i).setAreaCA(areaCA);
        }
    }

    
    
    private void setTrueAreaSize()
    {
        Area alpha;
        for(int i = 0; i < areas.size(); i++)
        {
            alpha = areas.get(i);
            double size = 0.0;
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int gridCellID = alpha.cells.get(j);
                int pos[] = pf.getPosition(gridCellID, dimX, dimY);
                double degreeLat = currentLatitudes[pos[0]];
                size+= Math.abs(Math.cos(Math.toRadians(degreeLat)));
            }
            areas.get(areas.indexOf(alpha)).trueSize = size;
        }        
    }

    
    //for a given latitude dimension, loads latitudeArray to current latitudes 
    //you can add more than the 2.0 degreee latitude resolution here.
    private void setCurrentLatitudes()
    {
        if(this.dimX == 90)
            currentLatitudes = latitudes90;
    }
    
    //resizes each value of the time series by the real size of the grid cell (cos(pji_i)).
    private void resizeTimeSeries()
    {
        for(int i = 0; i < dimY; i++)
        {
            for(int j = 0; j < dimX; j++)
            {
                double lat = currentLatitudes[j];
                double weight = Math.abs(Math.cos(Math.toRadians(lat)));
                for(int z = 0; z < dimT; z++)
                    dataArray[j][i][z] = dataArray[j][i][z]*weight;
            }
        }
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
            dataArray = new double[dimX][dimY][dimT];
            
            for(int i = 0; i < dimX; i++)
            {
                int timePos = 0, lonPos = 0;;
                for(int j = 0; j < temp[i].length; j++)
                {
                    dataArray[i][lonPos][timePos] = temp[i][j];
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
        }
        catch(IOException io)
        {
            System.out.println("Error parsing .mat file.");
        }
    }
    
    //Var. export functions go here
    public void setAreaCentroids()
    {
        double longitudes[] = new double[144];
        for(int i = 1; i < longitudes.length; i++)
            longitudes[i] = longitudes[i-1]+2.5;
        for(int i = 0; i < areas.size(); i++)
        {
            int[] center = getAreaCentroid(areas.get(i));
            double lat = latitudes90[center[0]];
            double lon = longitudes[center[1]];
            areas.get(i).cLat = lat;
            areas.get(i).cLon = lon;
        }
       
    }
    
    public int[] getAreaCentroid(Area area)
    {
        int cX,cY;
        double minDistance = Double.MAX_VALUE;
        int bestPos[] = new int[2];
        for(int i = 0; i < area.cells.size(); i++)
        {
            int cellID = area.cells.get(i);
            int posFrom[] = pf.getPosition(cellID, dimX, dimY);
            double avgDistance = 0.0;
            for(int j = 0; j < area.cells.size(); j++)
            {
                int cellID2 = area.cells.get(j);
                int posTo[] = pf.getPosition(cellID2, dimX, dimY);
                avgDistance+=euclideanDist(posFrom[0], posFrom[1], posTo[0], posTo[1]);
            }
            if(avgDistance <= minDistance)
            {
                minDistance = avgDistance;
                bestPos = posFrom;
            }
        }
        return bestPos;
    }
    
    public double euclideanDist(double x1, double y1, double x2, double y2)
    {
      return Math.sqrt(  Math.pow(x1-x2, 2)+Math.pow(y1-y2, 2));
    }
    
    private void printArea(Area area)
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"ensoArea.txt"));
            double[][] map = new double[dimX][dimY];
            for(int i = 0; i < area.cells.size(); i++)
            {
                int[] pos = pf.getPosition(area.cells.get(i), dimX, dimY);
                map[pos[0]][pos[1]] = 1;
            }
            for(int i = 0; i < dimX; i++)
            {
                for(int j = 0; j < dimY; j++)
                    out.print(map[i][j]+" ");
                out.println();
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void printAreaMap()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"map_areas.txt"));
            double[][] map = new double[dimX][dimY];
            for(int i = 0; i < areas.size(); i++)
            {
                Area alpha = areas.get(i);
                double color = alpha.color;
                for(int j = 0; j < alpha.cells.size(); j++)
                {
                    int pos[] = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                    if(map[pos[0]][pos[1]] > 0)//overlap
                        map[pos[0]][pos[1]] = 1.4;
                    else
                         map[pos[0]][pos[1]] = color;
                }
            }
            for(int i = 0; i < dimX; i++)
            {
                for(int j = 0; j < dimY; j++)                
                    out.print(map[i][j]+" ");
                out.println();
            }
            out.close();
            
        }
        catch(IOException ioe)
        {}
    }
    
    private void mapColoring(Vector<Area> areas)
    {
        Hashtable colorgraph[] = new Hashtable[areas.size()];
        Hashtable<Integer,Integer> mapIdCurrentToNew = new Hashtable<Integer,Integer>();
        Hashtable<Integer,Integer> mapIdNewToCurrent = new Hashtable<Integer,Integer>();
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            mapIdCurrentToNew.put(alpha.areaID, i);
            mapIdNewToCurrent.put(i, alpha.areaID);
        }
        for(int i = 0; i < colorgraph.length; i++)
            colorgraph[i] = new Hashtable<Integer,Integer>();
        //transform the map to a graph
        //a link will connect to areas in the map if they are adjacent
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            ArrayList<Integer> neighbors = findNeighborsMapColoring(areas,alpha);
            Integer from = mapIdCurrentToNew.get(alpha.areaID);
            for(int j = 0; j < neighbors.size(); j++)
            {
                Integer to = mapIdCurrentToNew.get(neighbors.get(j));
                colorgraph[from].put(to, to);
            }            
        }
        //I need it also as a Vector of areas...
        ArrayList<Area> areasMap = new ArrayList<>();
        for(int i = 0; i < colorgraph.length; i++ )
        {
            Area alpha = new Area(i);
            alpha.mapColEdges = new ArrayList<>();
            java.util.Enumeration<Integer> neighs = colorgraph[i].keys();
            while(neighs.hasMoreElements())
            {
                Integer neigh = neighs.nextElement();
                alpha.mapColEdges.add(new Edge(i, neigh));
            }
            areasMap.add(alpha);
        }
        ArrayList<Area> removeOrder = new ArrayList<>();
        AreaDegreeComparator degreeComp = new AreaDegreeComparator();
        while(!areasMap.isEmpty())
        {
            //sort the areas by degree in ascending order
            java.util.Collections.sort(areasMap,degreeComp);
            //remove the area with the lowest degree
            Area alpha = areasMap.remove(0);
            //check if it has degree less than five, if not we have a problem
            if(alpha.edges.size() > 4)
                System.out.println("Coloring problem");
            //remove all edges from alpha to other areas...
            for(int i = 0; i < alpha.edges.size(); i++)
            {
                Edge e = alpha.mapColEdges.get(i);
                for(int j = 0; j < areasMap.size(); j++)
                    if(areasMap.get(j).edges.contains(e))
                        areasMap.get(j).edges.remove(e);
            }
            removeOrder.add(alpha);
        }
        //initialize again the vector of areas, we will use it to set their colors
        areasMap = new ArrayList<>();
        for(int i = 0; i < colorgraph.length; i++ )
        {
            Area alpha = new Area(i);
            alpha.mapColEdges = new ArrayList<>();
            areasMap.add(alpha);
            //no need to know the edges here.        
        }
        double colorArray[] = {0.2,0.4,0.6,0.8,1.0};
        //now start removing the nodes in the remove order "stack-wise"
        while(!removeOrder.isEmpty())
        {
            Area alpha = removeOrder.remove(0);
            //and pick a color for alpha that is not used by its neighbors;
            boolean[] colorUsed = new boolean[5];
            for(int j = 0; j < colorUsed.length; j++)
                colorUsed[j] = false;
            java.util.Enumeration<Integer> neighs = colorgraph[alpha.areaID].keys();
            while(neighs.hasMoreElements())
            {
                Integer neigh = neighs.nextElement();
                Area dummy = new Area(neigh);
                dummy = areasMap.get(areasMap.indexOf(dummy));
                if(dummy.color > 0.0)
                {
                    if(dummy.color == 0.2)
                        colorUsed[0] = true;
                    else if(dummy.color == 0.4)
                        colorUsed[1] = true;
                    else if(dummy.color == 0.6)
                        colorUsed[2] = true;
                    else if(dummy.color == 0.8)
                        colorUsed[3] = true;
                    else colorUsed[4] = true;                                            
                }
            }
            //now scan the colorUsed array and pick a color for you...
            //one that is not used by your neighbors.
            //if you can't find such a color... problem             
            int colorid = -1;
            for(int j = 0; j < colorUsed.length; j++)
                if(!colorUsed[j])
                {
                    colorid = j;
               //     System.out.println(colorid);
                    break;
                }
            if(colorid == -1)
            {
                areasMap.get(areasMap.indexOf(alpha)).color = 1.2;
                System.out.println("Problem, no available colors.");
            }
            else
            {
                areasMap.get(areasMap.indexOf(alpha)).color = colorArray[colorid];
            }
                
        }
        //finally set the initial area colors
        for(int i = 0; i < areasMap.size(); i++)
        {
            Area alpha = areasMap.get(i);                        
            Integer oldId = mapIdNewToCurrent.get(alpha.areaID);
            Area dummy = new Area(oldId);
            
            this.areas.get(this.areas.indexOf(dummy)).color = alpha.color;
            
        }
    }
    
    private ArrayList<Integer> findNeighborsMapColoring(Vector<Area> areas,Area alpha)
    {
        
        ArrayList<Integer> neighIds = new ArrayList<>();
        for(int i = 0; i < alpha.cells.size(); i++)
        {
            Vector<Integer> cellNeighs = pf.getNeighborsCross(alpha.cells.get(i), dimX, dimY);
            for(int j = 0; j < cellNeighs.size(); j++)
            {
                Integer cellID = new Integer(cellNeighs.get(j));
                if(!alpha.cells.contains(cellID))
                {                    
                    for(int x = 0; x < areas.size(); x++)
                    {
                        Area beta = areas.get(x);
                        if(beta.cells.contains(cellID))
                        {                            
                            if(!neighIds.contains(beta.areaID))
                            {               
                                neighIds.add(new Integer(beta.areaID));
                            }
                        }
                    }
                }
            }
        }
        return neighIds;
    }

    private void exportForARI()
    {
        try
        {
            System.out.println("Exporting maps for ARI.");
          int map[][] = new int[dimX][dimY];
          int startID = 1;
          for(int i = 0; i < areas.size(); i++)
          {
              Area alpha = areas.get(i);
              for(int j = 0; j < alpha.cells.size(); j++)
              {
                  int pos[] = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                  map[pos[0]][pos[1]] = startID;
              }
              startID++;
          }
          PrintWriter out = new PrintWriter(new FileWriter(folderOut+"mapForARI.txt"));
          for(int i = 0; i < dimX; i++)
          {
              for(int j = 0; j < dimY; j++)
                  out.print(map[i][j]+" ");
              out.println();
          }
          out.close();
        }
        catch(IOException ioe)
        {}
    }


    private void exportLinkMap(Area from,String delimiter)
    {
        String linkMap[][] = new String[dimX][dimY];
        for(int i = 0; i < dimX; i++)
            for(int j = 0; j < dimY; j++)
                linkMap[i][j]=Double.toString(0.0);
        //out links
        for(int i = 0; i < from.edges.size(); i++)
        {
            AreaEdge e = from.edges.get(i);
            Area areaTo = new Area(e.to);
            areaTo = areas.get((areas.indexOf(areaTo)));
            for(int j = 0; j < areaTo.cells.size(); j++)
            {
                int areaToCellID = areaTo.cells.get(j);
                int pos[] = pf.getPosition(areaToCellID, dimX, dimY);
                linkMap[pos[0]][pos[1]] = Double.toString(e.weight);
            }
        }
        //in links
        for(int x = 0; x < areas.size(); x++)
        {
            Area area = areas.get(x);
            for(int i = 0; i < area.edges.size(); i++)
            {
                AreaEdge e = area.edges.get(i);
                if(e.to == from.areaID)
                {
                   for(int j = 0; j < area.cells.size(); j++) 
                   { 
                       int areaToCellID = area.cells.get(j);
                       int pos[] = pf.getPosition(areaToCellID, dimX, dimY);
                       linkMap[pos[0]][pos[1]] = Double.toString(e.weight);
                   }
                }
            }
        }
        //your position
        for(int i = 0; i < from.cells.size(); i++)
        {
            int cellidFrom = from.cells.get(i);
            int pos[] = pf.getPosition(cellidFrom, dimX, dimY);
            linkMap[pos[0]][pos[1]] = "NaN";
        }
        
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"linkMap_"+from.areaID+"_"+delimiter+".txt"));
            
            for(int i = 0; i < dimX; i++)
            {
                for(int j = 0; j < dimY; j++)
                    out.print(linkMap[i][j]+" ");
                out.println();
            }
            out.close();
            
            
        }
        catch(IOException ioe)
        {}

    }
    
    private void exportAreaIDAndCentroid()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"areaIDandCentroid.txt"));
            for(int i = 0; i < areas.size(); i++)
            {
                Area alpha = areas.get(i);
                out.println(alpha.areaID+"\t"+alpha.cLat+"\t"+alpha.cLon);
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }

    private void exportLagInfo(Area from, String delimiter)
    {
        
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"areaCentroidAndLag_"+delimiter+".txt"));
            HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
            for(int i = 0; i < areas.size(); i++)
                edges.addAll(areas.get(i).edges);
            Iterator<AreaEdge> it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                Area area = null;
                int lagSign = 0;
                if(e.to == from.areaID)//then this area points to you
                {
                    area = new Area(e.from);
                    area = areas.get(areas.indexOf(area));
                    lagSign = -1;
                }
                else if(e.from == from.areaID)//then you point to an area
                {
                    area = new Area(e.to);
                    area = areas.get(areas.indexOf(area));
                    lagSign = 1;
                }
                if(area!=null)    
                    out.println(area.cLat+" "+area.cLon+" "+(lagSign*Math.abs(e.lag)));
                
                
            }
            out.close();
        }
        catch(IOException ioe)
        {}                
    }

    public void exportGraphCLM(Vector<Area> areas, String fileDelim)
    {
        try
        {
             
            
            
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"finalGraph"+fileDelim+".clm"));
            out.println("Vertices\t"+areas.size());
            int numOfEdges = 0;
            Area alpha;
            for(int i = 0; i < areas.size(); i++)
            {
                alpha = areas.get(i);
                numOfEdges+=alpha.edges.size();
                out.println(alpha.areaID+"\t");
                for(int j = 0; j < alpha.cells.size(); j++)
                {
                    int cellid = alpha.cells.get(j);
                    out.print(cellid+" ");
                }
                out.print("\r\n");
            }//Finished with nodes.
            out.println("Edges\t"+numOfEdges);
            for(int i = 0; i < areas.size(); i++)
            {
                alpha = areas.get(i);
                for(int j = 0; j < alpha.edges.size(); j++)
                {
                    out.println(alpha.edges.get(j).from+"\t"+alpha.edges.get(j).to+"\t"+(alpha.edges.get(j).weight));
                }
            }
            out.close();
        }
        catch(IOException ioe)
        {
            ioe.printStackTrace();
        }
    }
    
    private void printStrengthMap()
    {
        double[][] map = new double[dimX][dimY];
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            double strength = alpha.strength;
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int pos[] = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                if(strength > map[pos[0]][pos[1]] )
                    map[pos[0]][pos[1]] = strength;
                //to plot the weaker area insted, comment the two lines above
                //and uncomment the four lines below
                
                //if(map[pos[0]][pos[1]] == 0)
                //  map[pos[0]][pos[1]] = strength; 
                //else if(strength < map[pos[0]][pos[1]] )
                //   map[pos[0]][pos[1]] = strength;                     
            }
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("strengthMap.txt"));
            for(int i = 0; i < dimX; i++)
            {
                for(int j = 0; j < dimY; j++)
                    out.print(map[i][j]+" ");
                out.println();
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void exportGraphNet(Vector<Area> areas, String delimiter)
    {
        try
        {
            Hashtable<Integer,Integer> newmap  = new Hashtable<Integer, Integer>();
            int cc = 1;
            for(int i = 0; i < areas.size(); i++)
            {
                newmap.put(areas.get(i).areaID, cc);
                cc++;
            }
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"finalGraph_"+delimiter+".net"));
            out.println("*Vertices "+areas.size());
            for(int i = 0; i < areas.size(); i++)
                out.println(newmap.get(areas.get(i).areaID));
            out.println("*arcs");
            for(int i = 0; i < areas.size(); i++)
            {
                Area alpha = areas.get(i);
                for(int j = 0; j < alpha.edges.size(); j++)
                {
                    AreaEdge e = alpha.edges.get(j);
                    if(e.from == alpha.areaID)
                        out.println(newmap.get(e.from)+" "+newmap.get(e.to)+" "+e.weight);
                }
            }
            out.close();
            
        }
        catch(IOException ioe)
        {}
    }

}


class AreaDegreeComparator implements Comparator<Area>{
    @Override
    public int compare(Area alpha, Area beta)
    {
        if(alpha.edges.size() > beta.edges.size())
            return 1;
        else if(alpha.edges.size() < beta.edges.size())
            return -1;
        else return 0;
    }
}


class ComparatorAreaTrueSizeBarlett implements Comparator<Area>
{
    public int compare(Area alpha, Area beta)
    {
        if(alpha.trueSize > beta.trueSize)
            return -1;
        else if(alpha.trueSize < beta.trueSize)
            return 1;
        else return 0;
    }
}
