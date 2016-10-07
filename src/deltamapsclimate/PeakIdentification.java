package deltamapsclimate;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class PeakIdentification {
    
    int dimX, dimY, dimT;
    Vector<GridCell> gridCells;
    Vector<GridCell> peaks;
    double[][] scoreMap;
    int neighSizeK;
  
    NewPublicFunctions pf = new NewPublicFunctions();
    PointDistanceComparator pointDistComp = new PointDistanceComparator();
  
    String folderOut;
    boolean[][] maskArray;
    double[][][] dataArray;
    double delta;//area identification threshold
    
    public PeakIdentification(DeltaMapsClimate parent)
    {
        System.out.println("--------- Peak Identification Starts --------------");
        this.dimT = parent.dimT;
        this.dimX = parent.dimX;
        this.dimY = parent.dimY;
        this.folderOut = parent.folderOut;
        this.neighSizeK = parent.neighSizeK;
        this.maskArray = parent.maskArray;
        this.dataArray = parent.dataArray;        
        this.delta = parent.delta;
        //step 1. For each grid cell construct it's local neighborhood.
        constructLocalNeighborhoods();       
        //step 2. For each grid cell calculate the local correlation field (or grid cell's score)
        inferScores();
        //step 3. Peak detection (inference of local maxima)
        inferPeaks();
    }
    
        
    private void inferPeaks()
    {
        System.out.println("Peak inference started.");
        peaks = new Vector<GridCell>();
        boolean[][] peakMap = new boolean[dimX][dimY];
        int nGridCells = gridCells.size();
        for(int i = 0; i < nGridCells; i++)
        {
            GridCell gc = gridCells.get(i);
            int myPos[] = pf.getPosition(gc.id,dimX,dimY);
            double myScore = scoreMap[myPos[0]][myPos[1]];
            if(myScore > delta)//criterion no 1.
            {
                Vector<Integer> neighs = pf.getNeighbors(gc.id, dimX, dimY);
                boolean peak = true;
                for(int j = 0; j < neighs.size(); j++)
                {
                    if(neighs.get(j) != gc.id )//remember we have our-self in the enighborhood
                    {
                        int[] neighPos = pf.getPosition(neighs.get(j), dimX, dimY);

                        if(scoreMap[neighPos[0]][neighPos[1]] == 0)//no data
                        {
                            peak = false;
                            break;
                        }
                        if(myScore < scoreMap[neighPos[0]][neighPos[1]] )// no peak, break
                        {
                            peak = false;
                            break;
                        }
                    }

                }
                peakMap[myPos[0]][myPos[1]]=peak;
                if(peak)
                    peaks.add(gridCells.get(i));

            }
                                    
        }
        
        //done, prepare and print a map;
        double[][] mapOut = scoreMap;
        for(int i = 0; i < dimX; i++)
        {
            for(int j = 0; j < dimY; j++)
                if(peakMap[i][j])
                    mapOut[i][j]=-1;
        }
        printMap(mapOut, "peakMap.txt");
        
        System.out.println("Found "+peaks.size()+" peaks");
    }
    

    
    private void inferScores()
    {
        System.out.println("Score inference started.");
        scoreMap = new double[dimX][dimY];
        for(int i = 0; i < gridCells.size(); i++)
        {
            GridCell gc = gridCells.get(i);
            double score = getScore(gc);
            gridCells.get(i).peakScore = score;
            int[] pos = pf.getPosition(gc.id, dimX, dimY);
            scoreMap[pos[0]][pos[1]] = score;
        }
    }
    
    private double getScore(GridCell gc)
    {
        int[] cellIDs = new int[gc.localNeigh.size()];        
        for(int i = 0; i < gc.localNeigh.size(); i++)
            cellIDs[i] = gc.localNeigh.get(i);
        
        double score = 0.0;
        double denom = 0.0;
        for(int i = 0; i < cellIDs.length; i++)
        {
            int posA[] = pf.getPosition(cellIDs[i], dimX, dimY);
            for(int j = (i+1); j < cellIDs.length; j++)
            {
                int posB[] = pf.getPosition(cellIDs[j], dimX, dimY);
                score+= pearsonCorrel(dataArray[posA[0]][posA[1]], dataArray[posB[0]][posB[1]]);
                denom++;
            }
        }
        
        return (score/denom);
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
    
    
    private void constructLocalNeighborhoods()
    {
        System.out.println("Constructing local neigborhoods");
        gridCells = new Vector<GridCell>();
        for(int x = 0; x < dimX; x++)
        {
            for(int y = 0; y < dimY; y++)
            {
                //first check if grid cell is not masked
                if(!maskArray[x][y])
                {
                    int gridCellID = pf.arrayToId(x, y, dimY);
                    GridCell gc = new GridCell(gridCellID);
                    gc.localNeigh.addAll(getLocalNeighborhoodDistance(gridCellID));
                    gridCells.add(gc);
                }
            }
        }
        System.out.println("Constructed the neighborhood for "+gridCells.size()+" grid cells.");
        //now let's to a cleaning of the data, some grid cells might not have a neighborhood size equal to neighSizze
        System.out.println("Cleaning data.");
        Vector<GridCell> toRemove = new Vector<GridCell>();
        for(int i = 0; i < gridCells.size(); i++)
        {
            if(gridCells.get(i).localNeigh.size()!=this.neighSizeK+1)//+1 because it contains it self 
                toRemove.add(gridCells.get(i));
        }
        gridCells.removeAll(toRemove);
        System.out.println("Grid cells with neighborhood: "+gridCells.size());
    }

    private Vector<Integer> getLocalNeighborhoodDistance(int gridCellID)
    {

        Vector<Integer> neighborhood = new Vector<Integer>();
        neighborhood.add(gridCellID);//Step 1. add yourself to the neighborhood.
        int posStart[] = pf.getPosition(gridCellID, dimX, dimY); //mark your start position
        
        //calculate the distance to all grid cells
        Vector<Point> kClosestPoints = new Vector<Point>();
        for(int i = 0; i < dimX; i++)
        {
            for(int j = 0; j < dimY; j++)
            {
                //get the e id of the grid cell
                int id = pf.arrayToId(i, j, dimY);
                //if it is not yourself
                if(id!=gridCellID)
                {
                    double distance = euclidianDistance(posStart[0], posStart[1], i, j);
                    Point p = new Point(id);
                    p.distance = distance;
                    kClosestPoints.add(p);
                }
            }
        }
        
        //sort by distance in ascending order
        java.util.Collections.sort(kClosestPoints,pointDistComp);
        //add the K closest points
        if(kClosestPoints.size() < this.neighSizeK)//do not bother, return the neighborhood
            return neighborhood;
        else
        {
            for(int i = 0; i < this.neighSizeK; i++)
            {
                int cellid = kClosestPoints.get(i).cellid;
                int pos[] = pf.getPosition(cellid, dimX, dimY);
                //if the kth closest point is masked do not add it
                if(!maskArray[pos[0]][pos[1]])
                    neighborhood.add(cellid);
            }
        }
        return neighborhood;
    }
    
    
    private double euclidianDistance(double x1, double y1, double x2, double y2)
    {
        return Math.sqrt(Math.pow(x1-x2, 2)+ Math.pow(y1-y2, 2));
    }
    
    private void printMap(double[][] map, String name)
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(this.folderOut+name));
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

}



class Point
{
    int cellid;
    double distance;
    
    public Point(int id)
    {
        this.cellid = id;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(o.getClass()!=this.getClass())
            return false;
        Point p = (Point)o;
        return this.cellid == p.cellid;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 23 * hash + this.cellid;
        return hash;
    }
}


class PointDistanceComparator implements Comparator<Point>
{
    @Override
    public int compare(Point n1, Point n2)
    {
        if(n1.distance > n2.distance)
            return 1;
        else if(n1.distance < n2.distance)
            return -1;
        else return 0;
                                                
    }
}