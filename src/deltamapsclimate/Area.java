package deltamapsclimate;

import java.util.ArrayList;
import java.util.Vector;
/**
 *
 * @author Ilias.Fountalis
 */
public class Area {
    
    int areaID; //in this implementation an area id is the adj. matrix cell id.
    Vector<Integer> cells;//the cell id in the area is the adjacency matrix id
    double score,size,trueSize;
    Vector<Integer> peaks;
    
    double[] areaCA;
    Vector<AreaEdge> edges = new Vector<AreaEdge>();
    double strength;
    
    double cLat,cLon,color;
    ArrayList<Edge> mapColEdges = new ArrayList<>();//just for map coloring
    
            
    public Area(int areaID)
    {
        this.areaID = areaID;        
    }
    
    
    public void setSize()
    {
        this.size = cells.size();
    }
    
    public void initPeaks(Vector<Integer> peaks)
    {
        this.peaks = new Vector<Integer>();
        for(int i = 0; i < peaks.size(); i++)
            if(this.cells.contains(peaks.get(i)))
                this.peaks.add(peaks.get(i));
    }
    
    public void setAreaCA(double[] ts)
    {
        this.areaCA = new double[ts.length];
        for(int x = 0; x < ts.length; x++)
            areaCA[x] = ts[x];
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 73 * hash + this.areaID;
        return hash;
    }

    @Override
    public boolean equals(Object o)
    {
        Area alpha = (Area)o;
        if(!(alpha instanceof Area))
            return false;
        else if(alpha.areaID == this.areaID)
            return true;
        else return false;
    }
    
}

//Edge class for MAP COLORING ONLY
class Edge
{
    //edge id's are adjacency matrix id's
    int from, to;

    
    public Edge(int from, int to)
    {
        this.from = from;
        this.to = to;        
    }
    

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 67 * hash + this.from;
        hash = 67 * hash + this.to;
        return hash;
    }
    
    @Override
    public boolean equals(Object o)
    {
        Edge e = (Edge)o;
        if(!(e instanceof Edge))
            return false;
        else if(e.from == this.from && e.to == this.to)
            return true;
        else if(e.from == this.to && e.to == this.from)
            return true;//graph is undirected...
        else return false;
    }
    
}
