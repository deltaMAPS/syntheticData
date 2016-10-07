package deltamapsclimate;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class AreaIdentification {
    
    double[][][] dataArray;
    boolean[][] maskArray;
    int dimX,dimY,dimT;
    double delta;
    int expandFactor = 5;
    String folderOut;
    int cc = 0;
    
    Vector<Area> areas;
    Vector<Integer> peaks;
    
    NewPublicFunctions pf = new NewPublicFunctions();
    
    GridCellScoreComp  gridScoreComp = new GridCellScoreComp();
    AreaScoreComp areaScoreComp = new AreaScoreComp();
    
    int mapID = 0;
    
    public AreaIdentification(DeltaMapsClimate parent, Vector<GridCell> peakIdentificationOut)
    {
        System.out.println("--------- Area Identification Starts --------------");
        this.dataArray = parent.dataArray;
        this.maskArray = parent.maskArray;
        this.dimT = parent.dimT;
        this.dimX = parent.dimX;
        this.dimY = parent.dimY;
        this.delta = parent.delta;
        this.folderOut = parent.folderOut;
        //Step 1. Construct the initial seed areas and also keep separately the peaks.
        loadData(peakIdentificationOut);
        //to speed things up, precalculate area scores
        for(int i = 0; i < areas.size(); i++)
            areas.get(i).score = getScore(areas.get(i).cells);
                
        //main loop. merge and expand until you can not do merging or expansion for any area.        
        
        while(true)
        {
            //System.out.println("Merging starts, current areas: "+areas.size());            
            boolean merged = mergeAreas();            
            //System.out.println("Merging stops, current areas: "+areas.size());
            //System.out.println("Expanding areas");
            boolean expanded = expandAreas();
            if(!merged && !expanded)
                break;
        }        

        for(int i = 0; i < areas.size(); i++)
            areas.get(i).initPeaks(peaks);
        
        printBorderMap(areas, "final");
        printOverlapMap(areas, "final");
        System.out.println("Final areas: "+areas.size());
        printAreaSize();
        
        exportAreaScoreMap();
        
        HashSet<Integer> cellsin = new HashSet<Integer>();
        for(int i = 0; i < areas.size(); i++)
            cellsin.addAll(areas.get(i).cells);
        double notmasked = 0;
        for(int i = 0; i < dimX; i++)
            for(int j = 0; j < dimY; j++)
            {
                if(!maskArray[i][j])
                    notmasked++;
            }
        exportGraphCLM();
        System.out.println("Cells in areas: "+((double)cellsin.size())/notmasked+"%");
                
        
    }
    
    //iteratively search for pairs of areas that can be merged, if more than one pair exists
    //merge the one whose union maximizes the average cross correlation. Terminate when no 
    //further merging is possible. Return true if you merged at least one pair.
    private boolean mergeAreas()
    {
        boolean merged = false;
        while(true) //merge as many areas as you can.
        {
            //(1) for each area get it's overlapping or adjacent neighbors.
            //(2) keep the best pair of areas to merge (best pair of areas == largest score OR max overlap)
            double bestScore = -1.0;
            int areaMergeAID = -1, areaMergeBID = -1;
            for(int i = 0; i < areas.size(); i++)
            {
                Area areaA = areas.get(i);
                Vector<Integer> geoConnAreas = getGeoConnectedAreas(areaA);
                for(int j = 0; j < geoConnAreas.size(); j++)
                {
                    int areaBID = geoConnAreas.get(j);
                    if(areaA.areaID > areaBID) //simple trick if areaA overlaps with areaB then the oposite will hold just get the score for one pair
                    {
                        Area areaB = new Area(areaBID);
                        areaB = areas.get(areas.indexOf(areaB));
                        double   score = getMergedScore(areaA, areaB);
                        if(score > bestScore)
                        {//then that is a better pair of areas to merge
                            areaMergeAID = areaA.areaID;
                            areaMergeBID = areaB.areaID;
                            bestScore = score;
                           }
                    }
                }            
            }//end now you should have the best area to merge
            if(bestScore > delta)
            {   //then I can merge the two areas.
                //remopve them from the areas.
                merged = true;
                Area areaToMergeA = new Area(areaMergeAID);
                areaToMergeA = areas.remove(areas.indexOf(areaToMergeA));
                Area areaToMergeB = new Area(areaMergeBID);
                areaToMergeB = areas.remove(areas.indexOf(areaToMergeB));
                //construct the new area
                Area mergedArea = new Area(areaToMergeA.areaID);//arbitrary id, dont care.
                mergedArea.cells = new Vector<Integer>();
                mergedArea.cells.addAll(getUnion(areaToMergeA, areaToMergeB));
                mergedArea.score = bestScore;
                mergedArea.size = mergedArea.cells.size();
                areas.add(mergedArea);
            }
            else//stop merging
                break;
        }//finished merging (end while)
        return merged;
    }
    
    
    
    private double getMergedScore(Area alpha, Area beta)
    {
        //(1) get the union of the grid cells.
        Vector<Integer> union = getUnion(alpha, beta);
        double score = 0.0;
        int nCells = union.size();
        for(int i = 0; i < nCells; i++)
        {
            int[] posFrom = pf.getPosition(union.get(i), dimX, dimY);
            for(int j = (i+1); j < nCells; j++)
            {
                int[] posTo = pf.getPosition(union.get(j), dimX, dimY);
                score+=pearsonCorrel(dataArray[posFrom[0]][posFrom[1]], dataArray[posTo[0]][posTo[1]]);
            }
        }
        double denom = nCells*(nCells-1)/2;
        return score/denom;
    }
    
    
    private Vector<Integer> getGeoConnectedAreas(Area anArea)
    {
        Vector<Integer> geoConnectedAreas = new Vector<Integer>();
        
        //step 1. construct a set that will have this area's grid cels as well as their
        //imediate border grid cells.
        HashSet<Integer> extendedAreaCells = new HashSet<Integer>();
        extendedAreaCells.addAll(anArea.cells);
        //and add the extended neighborhood here
        int nCells = anArea.cells.size();
        for(int i = 0; i < nCells; i++)
        {
            int cell = anArea.cells.get(i);
            Vector<Integer> neighs = pf.getNeighborsCross(cell, dimX, dimY);
            extendedAreaCells.addAll(neighs);
        }
        //now search for overlaps
        int nAreas = areas.size();
        for(int i = 0; i < nAreas; i++)
        {
            Area alpha = areas.get(i);
            if(alpha.areaID!=anArea.areaID)
            {
                HashSet<Integer> union = new HashSet<Integer>();
                union.addAll(extendedAreaCells);
                union.addAll(alpha.cells);
                if(union.size() < (alpha.cells.size()+extendedAreaCells.size()))
                    geoConnectedAreas.add(alpha.areaID);
            }
        }
        return geoConnectedAreas;
    }
    
    private Vector<Integer> getUnion(Area alpha, Area beta)
    {
        HashSet<Integer> union = new HashSet<Integer>();
        union.addAll(alpha.cells);
        union.addAll(beta.cells);
        Vector<Integer> toReturn = new Vector<Integer>();
        toReturn.addAll(union);
        return toReturn;
    }

    
    /*
    Expand the area that (a) can be expanded and (b) has the highest score.
    if an area has been expanded check if you can merge. if yes return to merging loop
    else repeat expansion
    */
    private boolean expandAreas()
    {        
        
        boolean startMerging = false;
        boolean expanded = false;
        while(!startMerging)
        {
            expanded = false;
            for(int i = 0; i < areas.size(); i++)
                areas.get(i).score = getScore(areas.get(i).cells);
            //sort the areas by score in sescending order.
            java.util.Collections.sort(areas,areaScoreComp);
            //try to expand areas, starting from the area with the highest score                
            int nAreas = areas.size(); 
            Area expandedArea = null;
            for(int i = 0; i < nAreas; i++)
            {
                Area area = areas.get(i);
                //System.out.println("Expanding area: "+area.areaID+" score: "+area.score);
                expandedArea = expandArea(area);                
                if(expandedArea.cells.size() > area.cells.size() ) //then the area has been expanded
                {
                    //System.out.println("Expanded area with score: "+expandedArea.score);
                    //System.out.println("Expanded area size: "+expandedArea.cells.size());
                    //remove the old area and add the new one
                    expanded = true;
                    expandedArea.size = expandedArea.cells.size();
                    int index = areas.indexOf(area);
                    areas.remove(index);
                    areas.add(index, expandedArea);                    
                    Vector<Integer> neighIDs = getGeoConnectedAreas(expandedArea);
                    //System.out.println("Start merging?");
                    for(int j = 0; j < neighIDs.size(); j++)
                    {
                        Area foo = new Area(neighIDs.get(j));
                        foo = areas.get(areas.indexOf(foo));
                        double score = getMergedScore(expandedArea, foo);
                        if(score>delta)
                        {
                            startMerging=true;
                            //System.out.println("Ya.");
                            break;   
                        }
                        //else
                            //System.out.println("Nope");
                    }
                }          
                if(startMerging)
                    break;
            }
            if(!expanded)
                break;
        }
        return expanded;
    }
    

    
    
    private Area expandArea(Area area)
    {
        //initialize a set that will hold all grid cells in this area
        HashSet<Integer> inArea = new HashSet<Integer>(); 
        //and add all that are in already
        inArea.addAll(area.cells);
        //initialize the area
        Area expandedArea = new Area(area.areaID);
        
        
        expandedArea.cells = new Vector<Integer>();        
        //add the ones that are inside already
        expandedArea.cells.addAll(area.cells);
        expandedArea.score = area.score;//getScore(expandedArea.cells);
        //find the geographically connected neighbors
        Vector<GridCell> neighboringCells = new Vector<GridCell>(); //neighboring cell id's are adj. matrix id's
        for(int i = 0; i < expandedArea.cells.size(); i++)
        {
            int cellID = expandedArea.cells.get(i);
            Vector<Integer> tempNeighs = pf.getNeighborsCross(cellID, dimX, dimY);
            for(int j = 0; j < tempNeighs.size(); j++)
            {
                int neighID = tempNeighs.get(j);
                int pos[] = pf.getPosition(neighID, dimX, dimY);
                //if the neighbor is not masked and not already in the area
                if(maskArray[pos[0]][pos[1]] == false && !inArea.contains(neighID))
                {
                    GridCell n = new GridCell(neighID);
                    if(!neighboringCells.contains(n))
                        neighboringCells.add(n);
                }
            }
        }
        
        int nodesEntered = 0;
        //ready to start expanding area
        //Stop criterion I:  best node score is less than the threshold
        //Stop criterion II: when we do not have any neighboring cells left
        //Stop criterion III: when you've added enough nodes.
        while(neighboringCells.size() > 0) //Stop criterion II check here.
        {
            if(nodesEntered == expandFactor)
                break;//added enough nodes.
            //calculate the score of all the grid cells.
            for(int i = 0; i < neighboringCells.size(); i++)            
                if(neighboringCells.get(i).updateScore == false)//check if I already have the current node's score
                    neighboringCells.get(i).score = getScore(neighboringCells.get(i).id,expandedArea);            
            //sort the nodes by score
            java.util.Collections.sort(neighboringCells,gridScoreComp);
            //and get the node with the highest score
            GridCell bestCandidate = neighboringCells.remove(0);
                                    
            if(bestCandidate.score > delta)//Stop criterion I check here
            {
                double areaPairs = expandedArea.cells.size()*(expandedArea.cells.size()-1)/2;
                double currentAreaScore = areaPairs*expandedArea.score;
                double bcScore = bestCandidate.score*expandedArea.cells.size();
                expandedArea.score = (bcScore+currentAreaScore)/ (double)(expandedArea.cells.size()*(expandedArea.cells.size()+1)/2);
                
                //ok he passed add him to the area
                expandedArea.cells.add(bestCandidate.id);
                //make him unavailable 
                inArea.add(bestCandidate.id);
                //upadte the number of nodes in area
                nodesEntered++;
                //update the score of nodes
                for(int x = 0; x < neighboringCells.size(); x++)                                    
                {
                    GridCell toUpdate = neighboringCells.get(x);
                    neighboringCells.get(x).score = updateScore(toUpdate.id,bestCandidate.id, expandedArea.cells.size()-1,toUpdate.score);            
                    neighboringCells.get(x).updateScore=true;
                }
                                    
                //get his geographically connected neighbors and add them as candidates for expansion
                Vector<Integer> tempNeighs = pf.getNeighborsCross(bestCandidate.id, dimX, dimY);
                for(int i = 0; i < tempNeighs.size(); i++)
                {
                    int[] tempPos = pf.getPosition(tempNeighs.get(i), dimX, dimY);
                    if(maskArray[tempPos[0]][tempPos[1]] == false && !inArea.contains(tempNeighs.get(i)))
                    {
                        GridCell n = new GridCell(tempNeighs.get(i));
                        if(!neighboringCells.contains(n))
                                neighboringCells.add(n);
                    }

                }//finished adding new neighbors...
            }
            else break;
        }
        return expandedArea;//return the area
    }
    

    /*
    neighbor : the node that needs to update the score
    bestCandidate: the last node to enter the area
    areaSize: areaSize before the bestCandidate joined the area
    currentScore: currentScore of neighbor to all other cells in the area
    */
    private double updateScore(int neighbor, int bestCandidate,  int areaSize, double currentScore)
    {
        int [] posFrom = pf.getPosition(neighbor, dimX, dimY);
        int [] posTo = pf.getPosition(bestCandidate, dimX, dimY);
        double corr = pearsonCorrel(dataArray[posFrom[0]][posFrom[1]], dataArray[posTo[0]][posTo[1]]);
        return (areaSize*currentScore +corr)/(double)(areaSize+1);
    }
    
    private double getScore(int from ,Area areaTo)
    {
        double score = 0;
        int[] posFrom = pf.getPosition(from, dimX, dimY);        
        double[] tsFrom = dataArray[posFrom[0]][posFrom[1]];
        for(int i = 0; i < areaTo.cells.size(); i++)
        {
            int posTo[] = pf.getPosition(areaTo.cells.get(i), dimX, dimY);
            double[] tsTo = dataArray[posTo[0]][posTo[1]];
            score += pearsonCorrel(tsFrom, tsTo);
        }
        return score/=(double)areaTo.cells.size();
    }

    
    private double getScore(Vector<Integer> areaCells)
    {
        double score = 0;
        double denom = 0;
        int nCells = areaCells.size();
        for(int i = 0; i < nCells; i++)
        {
            int posFrom[] = pf.getPosition(areaCells.get(i), dimX, dimY);
            double[] tsFrom = dataArray[posFrom[0]][posFrom[1]];
            for(int j = (i+1); j < nCells; j++)
            {
                int posTo[] = pf.getPosition(areaCells.get(j), dimX, dimY);
                double[] tsTo = dataArray[posTo[0]][posTo[1]];
                score+=pearsonCorrel(tsFrom, tsTo);
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

    private void exportGraphCLM()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"finalAreas.clm"));
            out.println("Vertices\t"+areas.size());
            for(int i = 0; i < areas.size(); i++)
            {
                Area alpha = areas.get(i);
                out.println(alpha.areaID);
                for(int j = 0; j < alpha.cells.size(); j++)
                    out.print(alpha.cells.get(j)+" ");
                out.print("\r\n");
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
 
    
    
    private void loadData(Vector<GridCell> peakIdentificationOut)
    {
        areas = new Vector<Area>();
        peaks = new Vector<Integer>();
        int nSeedAreas = peakIdentificationOut.size();
        int areaID = 0;
        for(int i = 0; i < nSeedAreas; i++)
        {
            GridCell gc = peakIdentificationOut.get(i);
            peaks.add(gc.id);
            Area alpha = new Area(areaID); areaID++;
            alpha.cells = new Vector<Integer>();
            alpha.cells.addAll(gc.localNeigh);
            areas.add(alpha);
        }
        System.out.println("Initial seed areas: "+areas.size());
    }
    
    
    private void printAreaSize()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"areaSize.txt"));
            for(int i = 0; i < areas.size(); i++)
                out.println(areas.get(i).cells.size());
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void printOverlapMap(Vector<Area> areas, String delimeter)
    {
        int[][] map = new int[dimX][dimY];
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int pos[] = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                map[pos[0]][pos[1]]+=1;
            }
        }
        for(int i = 0; i < peaks.size(); i++)
        {
            int pos[] = pf.getPosition(peaks.get(i), dimX, dimY);
            map[pos[0]][pos[1]] = -1;
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"overlapMap_"+delimeter+".txt"));
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
    
    private void printOverlapMap2(Vector<Area> areas, String delimeter)
    {
        double[][] map = new double[dimX][dimY];
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int pos[] = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                if(map[pos[0]][pos[1]] == 0)
                    map[pos[0]][pos[1]]=alpha.areaID;
                else
                    map[pos[0]][pos[1]] = ((double)(map[pos[0]][pos[1]]+alpha.areaID))/2.0;
            }
        }
        for(int i = 0; i < peaks.size(); i++)
        {
            int pos[] = pf.getPosition(peaks.get(i), dimX, dimY);
            map[pos[0]][pos[1]] = -1;
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"overlapMap_"+delimeter+".txt"));
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

    
    private void printBorderMap(Vector<Area> areas, String delimeter)
    {
        int map[][] = new int[dimX][dimY];
        int borderID = (int)Math.floor(Math.random()*areas.size()+1);
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int cid = alpha.cells.get(j);
                Vector<Integer> neighPos = pf.getNeighborsCross(cid, dimX, dimY);
                int isBorder = 0; //yes if <4.
                for(int z = 0; z < neighPos.size(); z++)
                {
                    if(alpha.cells.contains(neighPos.get(z)))
                        isBorder++;
                }
                if(isBorder!=4)
                {
                    int pos[] = pf.getPosition(cid, dimX, dimY);
                    map[pos[0]][pos[1]]=borderID;
                }
            }
            borderID = (int)Math.floor(Math.random()*areas.size())+1;
        }
        for(int i = 0; i < peaks.size(); i++)
        {
            int pos[] = pf.getPosition(peaks.get(i), dimX, dimY);
            map[pos[0]][pos[1]] = -1;
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"borderMap_"+delimeter+".txt"));
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
    
    
    
    private void exportAreaScoreMap()
    {
        double[][] map = new double[dimX][dimY];
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            double score = getScore(alpha.cells);
            for(int j = 0; j < alpha.cells.size(); j++)
            {
                int[] pos = pf.getPosition(alpha.cells.get(j), dimX, dimY);
                map[pos[0]][pos[1]]=score;
                
            }
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(folderOut+"areaScoreMap.txt"));
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


class GridCellScoreComp implements Comparator<GridCell>
{
    @Override
    public int compare(GridCell n1, GridCell n2)
    {
        if(n1.score < n2.score)
            return 1;
        else if(n1.score > n2.score)
            return -1;
        else return 0;
                                                
    }
}


class AreaScoreComp implements Comparator<Area>
{
    @Override
    public int compare(Area n1, Area n2)
    {
        if(n1.score < n2.score)
            return 1;
        else if(n1.score > n2.score)
            return -1;
        else return 0;
                                                
    }
}


