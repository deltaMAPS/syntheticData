package deltamapsclimate;

import java.util.Vector;

/**
 *
 * @author ilias.fountalis
 */


public class GridCell {
        
        int id;
        double peakScore;
        Vector<Integer> localNeigh; //local neighborhood for peak identification
        double score; //avg. correlation to all other cells in the area that you're trying to enter.
        boolean updateScore = false;
        
       
        
        public GridCell(int id)
        {
            this.id = id;
            localNeigh = new Vector<Integer>();
        }
        
        @Override
        public boolean equals(Object o)
        {
            if(this.getClass() != o.getClass())
                return false;
            GridCell gc = (GridCell)o;
            return this.id == gc.id;
        }
                       
        @Override
        public int hashCode() {
            int hash = 5;
            hash = 67 * hash + this.id;
            return hash;
        }
}
