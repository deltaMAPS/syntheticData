package deltamapsclimate;


import java.util.Vector;
/**
 *
 * @author ilias.fountalis
 */
public class NewPublicFunctions {
    
    public final double EARTHRAD = 6371.009;//in km


    public NewPublicFunctions()
    {
    }
    
    public Vector<Integer> getNeighbors(int node, int dimX,int dimY)
    {
        node = node-1; //when the area nodes start counting from 1 to N, remove it if count is 0 to N-1

        Vector<Integer> neighbors = new Vector<Integer>();
        int []pos = getPosition(node+1,dimX,dimY);//cause in reality we start counting from 1 not 0.

        if(pos[0]+1 < dimX )
            neighbors.add(new Integer(arrayToId(pos[0]+1,pos[1],dimY)));
        if(pos[0]-1 >= 0 )
            neighbors.add(new Integer(arrayToId(pos[0]-1,pos[1],dimY)));
        if(pos[1]+1 < dimY )
            neighbors.add(new Integer(arrayToId(pos[0],pos[1]+1,dimY)));
        if(pos[1]-1 >= 0 )
            neighbors.add(new Integer(arrayToId(pos[0],pos[1]-1,dimY)));

        if(pos[1]+1 < dimY && pos[0]+1 < dimX)
            neighbors.add(new Integer(arrayToId(pos[0]+1,pos[1]+1,dimY)));
        if(pos[1]+1 < dimY && pos[0]-1 >= 0)
            neighbors.add(new Integer(arrayToId(pos[0]-1,pos[1]+1,dimY)));
        if(pos[1]-1 >= 0 && pos[0]+1 < dimX)
            neighbors.add(new Integer(arrayToId(pos[0]+1,pos[1]-1,dimY)));
        if(pos[1]-1 >= 0 && pos[0]-1 >= 0)
            neighbors.add(new Integer(arrayToId(pos[0]-1,pos[1]-1,dimY)));

        if(pos[1]+1 == dimY)
            neighbors.add(new Integer(arrayToId(pos[0], 0,dimY)));
        if(pos[1]-1 < 0)
            neighbors.add(new Integer(arrayToId(pos[0],dimY-1,dimY)));
        if(pos[1]+1 == dimY && pos[0]-1 >= 0)
            neighbors.add(new Integer(arrayToId(pos[0]-1, 0,dimY)));
        if(pos[1]+1 == dimY && pos[0]+1 < dimX)
            neighbors.add(new Integer(arrayToId(pos[0]+1, 0,dimY)));
        if(pos[1]-1 < 0 && pos[0]+1 < dimX)
            neighbors.add(new Integer(arrayToId(pos[0]+1, dimY-1,dimY)));
        if(pos[1]-1 < 0 && pos[0]-1 >= 0)
            neighbors.add(new Integer(arrayToId(pos[0]-1, dimY-1,dimY)));

        return neighbors;
    }


    //position of a node in an array can be a short value...
    public int[] getPosition(int id,int dimX, int dimY)
    {
        int[] pos =  new int[2];
        if(id%dimY == 0)
        {
            pos[0] = (id/dimY-1);
            pos[1] = (dimY-1);
        }
        else
        {
            pos[0] = (int)Math.floor(id/dimY);
            pos[1] = (id%dimY-1);
        }
        return pos;
    }
    
    
    //neighbor ids need to be of Integer value
    public Vector<Integer> getNeighborsCross(int node, int dimX,int dimY)
    {
        Vector<Integer> neighbors = new Vector<Integer>();        
        int []pos = getPosition( node, dimX, dimY);
        
        if(pos[0]+1 < dimX )
            neighbors.add(new Integer(arrayToId((pos[0]+1),pos[1],dimY)));
        if(pos[0]-1 >= 0 )
            neighbors.add(new Integer(arrayToId((pos[0]-1),pos[1],dimY)));
        if(pos[1]+1 < dimY )
            neighbors.add(new Integer(arrayToId(pos[0],(pos[1]+1),dimY)));
        if(pos[1]-1 >= 0 )
            neighbors.add(new Integer(arrayToId(pos[0],(pos[1]-1),dimY)));

        
        if(pos[1]+1 == dimY)
            neighbors.add(new Integer(arrayToId(pos[0], 0,dimY)));
        if(pos[1]-1 < 0)
            neighbors.add(new Integer(arrayToId(pos[0],(dimY-1),dimY)));        
        
        return neighbors;
    }

    public int arrayToId(int x, int y,int dimY)
    {
        return (dimY*x +(y+1));
    }
    
    public double distance3(double lat_a, double lon_a,double lat_b, double lon_b)
    {
        double  dLat = lat_a-lat_b;
        double dLon = lon_a-lon_b;
        double alpha = Math.sin(dLat/2) * Math.sin(dLat/2) + Math.cos(lat_a)*Math.cos(lat_b)*Math.sin(dLon/2)*Math.sin(dLon/2);
        double c = 2 * Math.atan2(Math.sqrt(alpha), Math.sqrt(1-alpha));
        return EARTHRAD*c;

        //return EARTHRAD*Math.acos(Math.sin(a.y)*Math.sin(b.y)+ Math.cos(a.y)*Math.cos(b.y)*Math.cos((a.x-b.x)));

    }


}

