package deltamapsclimate;

/**
 *
 * @author ilias.fountalis
 */
public class AreaEdge {
    
    int from, to; 
    double weight;
    double sij;
    
    //double minLag, maxLag;
    double lag;
    boolean undirected = false;
    double[] correlogram;
    double[] sigcorrs;
    double[] bvar;
    int minLag,maxLag;
    
    public AreaEdge()
    {}
    
    public AreaEdge(int from, int to, double weight)
    {
        this.from = from;
        this.to = to;
        this.weight = weight;
    }
    
    

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 11 * hash + this.from;
        hash = 11 * hash + this.to;
        return hash;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(this.getClass() != o.getClass())
            return false;
        AreaEdge e = (AreaEdge)o;
        return this.from == e.from && this.to == e.to;
    }
}
