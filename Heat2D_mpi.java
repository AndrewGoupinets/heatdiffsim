import java.util.Date;
import mpi.*;

public class Heat2D_mpi
{
    private static double a = 1.0;  // heat speed
    private static double dt = 1.0; // time quantum
    private static double dd = 2.0; // change in system
    private static int rank;        // rank of this machine
    private static int startPos;    // starting position of square for machine
    private static int endPos;      // ending position of square for machine
    
    //used for message passing
    final static int id = 1;
    



    /*
        main
        Starting point of the program. Initializes variables and goes through  
        the logic of the program
    */
    public static void main( String[] args ) throws MPIException
    {
        // Start the MPI library.
        MPI.Init( args );
        
	// verify arguments
	if ( args.length != 4 )
        {
	    System.out.println( "usage: " + 
			 "java Heat2D size max_time heat_time interval" );
	    System.exit( -1 );
	}

        
        //initializes with original 4 arguments
	int size = Integer.parseInt( args[0] );
	int max_time = Integer.parseInt( args[1] );
	int heat_time = Integer.parseInt( args[2] );
	int interval  = Integer.parseInt( args[3] );
	double r = a * dt / ( dd * dd );
        
        
        //determine rank
        rank = MPI.COMM_WORLD.Rank();
        
	// create a space
	double[][][] z = new double[2][size][size];
	for ( int p = 0; p < 2; p++ )
        {
            for ( int x = 0; x < size; x++ )
            {
        	for ( int y = 0; y < size; y++ )
                {
                    z[p][x][y] = 0.0; // no heat or cold
                }    
            }
        }

	
        int slice = size / MPI.COMM_WORLD.Size();
        int remainder = size % MPI.COMM_WORLD.Size();
        
        startPos = slice * rank;
        endPos = startPos + slice - 1;
        
        if(rank == MPI.COMM_WORLD.Size() - 1 && remainder != 0)
        {
            endPos = size - 1;
        }
        
	// start a timer
	Date startTime = new Date( );
	
	// simulate heat diffusion
	for ( int t = 0; t < max_time; t++ )
        {
	    int p = t % 2; // p = 0 or 1: indicates the phase
	    
            
            //the first three loops to be done on all machines
	    universalLoops(z, p, size, t, heat_time);
            

            //share boundary information between machines
            shareBoundaryData(size, z, p);
            

	    // updates the Master z[][][] with all other's
            //machine's slices
            updateMaster(size, z, p);


	    // display intermediate results
            if(rank == 0)
            {
                printMaster(interval, t, max_time, size, z, p);
            }
            
            
	    // perform forward Euler method
	    forwardEulerMethod(p, size, z, r);
	} // end of simulation
	
	// finish the timer
            if(rank == 0)
            {
                Date endTime = new Date( );
                System.out.println( "Elapsed time = " + 
		( endTime.getTime( ) - startTime.getTime( ) ) );
            }
        
	// Terminate the MPI library.
        MPI.Finalize( );
    }
    



    /*
        forwardEulerMethod
        Computes the next state of z[][][] by using the neighboring elements
        to each element and saves them on the non-active dimension 
        (if p=1, p2 = 0, and vice-versa)
    */
    public static void forwardEulerMethod(int p, int size, double[][][] z,
            double r)
    {
        int p2 = (p + 1) % 2;
        
        try
        {
            //if there is only one machine, this machine is both first and last
            if(MPI.COMM_WORLD.Size() == 1)
            {
                for ( int h = 1; h < size - 1; h++ )
                {
                    for ( int v = 1; v < size - 1; v++ )
                    {
                        z[p2][h][v] = z[p][h][v] + 
                        r * ( z[p][h+1][v] - 2 * z[p][h][v] + z[p][h-1][v] ) +
                        r * ( z[p][h][v+1] - 2 * z[p][h][v] + z[p][h][v-1] );
                    }  
                }
            }
            //if first machine, you have to ignore the first element
            else if(rank == 0)
            {
                for ( int h = 1; h <= endPos; h++ )
                {
                    for ( int v = 1; v < size - 1; v++ )
                    {
                        z[p2][h][v] = z[p][h][v] + 
                        r * ( z[p][h+1][v] - 2 * z[p][h][v] + z[p][h-1][v] ) +
                         r * ( z[p][h][v+1] - 2 * z[p][h][v] + z[p][h][v-1] );
                    }  
                }
            }
            //if last machine, you have to ignore the last element
            else if(rank == MPI.COMM_WORLD.Size() - 1)
            {
                for ( int h = startPos; h < size - 1; h++ )
                {
                    for ( int v = 1; v < size - 1; v++ )
                    {
                        z[p2][h][v] = z[p][h][v] + 
                        r * ( z[p][h+1][v] - 2 * z[p][h][v] + z[p][h-1][v] ) +
                        r * ( z[p][h][v+1] - 2 * z[p][h][v] + z[p][h][v-1] );
                    }  
                }
            }
            //else you are a middle machine and don't need to worry
            else
            {
                for ( int h = startPos; h <= endPos; h++ )
                {
                    for ( int v = 1; v < size - 1; v++ )
                    {
                        z[p2][h][v] = z[p][h][v] + 
                        r * ( z[p][h+1][v] - 2 * z[p][h][v] + z[p][h-1][v] ) +
                        r * ( z[p][h][v+1] - 2 * z[p][h][v] + z[p][h][v-1] );
                    }  
                }
            }
        }
        catch(MPIException e)
        {
            System.out.println("An error occured in forwardEulerMethod");
        }
    }
    



    /*
        printMaster
        Prints out z[][][] after every interval or right before max_time 
        is reached. Runs a small math function to divide every number by
        2 and floor them as ints
    */
    public static void printMaster(int interval, int t, int max_time, int size, 
            double[][][] z, int p)
    {
        if(interval != 0 && 
		 ( t % interval == 0 || t == max_time - 1 ) )
        {
            System.out.println( "time = " + t );
            for ( int y = 0; y < size; y++ )
            {
                for ( int x = 0; x < size; x++ )
                {
                    System.out.print( (int)( Math.floor(z[p][x][y] / 2) ) 
					  + " " );
                }
            System.out.println( );
            }
	System.out.println( );
	}
    }
    



    /*
        updateMaster
        Updates the master machine's z[][][] double array by receiving every
        other machine's slice
        Master machine communicates with the other machine's in order, inquiring
        about their start and end positions, generating a buffer based on that,
        and filling that buffer with their slices
    */
    public static void updateMaster(int size, double[][][] z, int p)
    {
        try
        {
            //if machine is the master, inquire for all other machine's slices
            if(rank == 0)
            {
                //Master inquires for slice information from other machines
                for(int i = 1; i < MPI.COMM_WORLD.Size(); i++)
                {
                    //asks machine rank i for it's start and end positions
                    int[] startEnd = new int[2];
                    MPI.COMM_WORLD.Send(startEnd, 0, 2, MPI.INT, i, id);

                    //receives machine rank i's start and end positions
                    MPI.COMM_WORLD.Recv(startEnd, 0, 2, MPI.INT, i, id);
                    
                    //uses this information to generate a buffer of an appropriate
                    //size
                    int bufferSize = (startEnd[1] - startEnd[0] + 1) * size;

                    //System.out.println("bufferSize: " + bufferSize);

                    double[] buffer = new double[bufferSize];

                    //asks machine rank i for it's slice information
                    MPI.COMM_WORLD.Send(buffer, 0, bufferSize, MPI.DOUBLE, i, id);

                    //receive's machine rank i's slice information
                    MPI.COMM_WORLD.Recv(buffer, 0, bufferSize, MPI.DOUBLE, i, id);
                          
                    //modifies master's z[][][] with buffer
                    int counter = 0;
                    for(int h = startEnd[0]; h <= startEnd[1]; h++)
                    {
                        for(int v = 0; v < size; v++)
                        {
                            z[p][h][v] = buffer[counter];
                            counter++;
                        }
                    }
                }
            }
            
            //all other machines send information to master
            else
            {
                //generates int[] with start and end positions
                int[] startEnd = new int[]{startPos, endPos};
                int[] temp = new int[2];
                
                
                //waits for master to ask for start and end positions
                MPI.COMM_WORLD.Recv(temp, 0, 2, MPI.INT, 0, id);
                
                //sends start and end positions to master
                MPI.COMM_WORLD.Send(startEnd, 0, 2, MPI.INT, 0, id);
                
                //generates array representing it's values
                int bufferSize = (endPos - startPos + 1) * size;
                double[] buffer = new double[bufferSize];
                
                //waits to be asked for their slice information
                MPI.COMM_WORLD.Recv(buffer, 0, bufferSize, MPI.DOUBLE, 0, id);
                
                //populates said array
                int count = 0;
                for(int h = startPos; h <= endPos; h++)
                {
                    for(int v = 0; v < size; v++)
                    {
                        buffer[count] = z[p][h][v];
                        count++;
                    }
                }
                
                //sends slice information
                MPI.COMM_WORLD.Send(buffer, 0, bufferSize, MPI.DOUBLE, 0, id);
            }  
        }
        catch(MPIException e)
        {
            System.out.println("An error occured in updateMaster");
        }    
    }




    /*
        universalLoops
        These loops are to be run in all the machines as applicable.
        These simulate the heat be diffused and consists of three functions.
    */
    public static void universalLoops(double[][][] z, int p, int size, int t,
            int heat_time)
    {
        //FIRST LOOP
        // two left-most and two right-most columns are identical
	// this loop is to be handled with the outermost machines, or by
        // only the master machine if only one machine is in use
        for ( int y = 0; y < size; y++ )
        {
            try
            {
                //check for master first
                if(rank == 0)
                {
                    z[p][0][y] = z[p][1][y];
                }
                //check for last machine
                //could also be the master if only one machine
                if(rank == MPI.COMM_WORLD.Size() - 1)
                {
                    z[p][size - 1][y] = z[p][size - 2][y];
                }
            }
            catch(MPIException e)
            {
                System.out.println("An error occurred in the first loop");
            }
	}


        //SECOND LOOP
	// two upper and lower rows are identical
        // every machine goes through it's slice and changes the top two
        // and bottom two rows
        for(int i = startPos; i <= endPos; i++)
        {
            z[p][i][0] = z[p][i][1];
            z[p][i][size - 1] = z[p][i][size -2];
        }


        //THIRD LOOP
	// keep heating the bottom until t < heat_time
        try
        {
            if(t < heat_time)
            {
                //odd # of computers
                if(MPI.COMM_WORLD.Size() % 2 == 1)
                {
                     //if machine is the middle most one
                     if(rank == MPI.COMM_WORLD.Size() / 2)
                     {
                         for ( int x = size /3; x < size / 3 * 2; x++ )
                         {
                            z[p][x][0] = 19.0; // heat
                         }

                     }
                 }
             //even # of computers
             else
             {
                //selects two innermost machines
                if(rank == MPI.COMM_WORLD.Size() / 2)
                {
                    if ( t < heat_time )
                    {
                        for ( int x = startPos; x < size / 3 * 2; x++ )
                        {
                            z[p][x][0] = 19.0; // heat
                        }
                    }
                    }
                    else if(rank == MPI.COMM_WORLD.Size() / 2 - 1)
                    {
                        for ( int x = size / 3; x <= endPos; x++ )
                        {
                            z[p][x][0] = 19.0; // heat
                        }
                    }
                }
            }
        }
        catch(MPIException e)
        {
            System.out.println("An error occured in the third loop");
        }
    }




    /*
        shareBoundaryData
        Every machine sends and sets the boundary data immediately outside of 
        it's start and end positions as applicable. This is accomplished by 
        dividing the machines into two groups, one that first sends and then 
        receives, the other that first receives then sends
    */
    public static void shareBoundaryData(int size, double[][][] z, int p)
    {
            //master machines sends boundary data to the right and receives some
            //from the right
        try
        {
            double[] workingColumn = new double[size]; // double[] used to send
                                                       // and receive boundaries
            
            
            //if there is only one machine working, no need to share data
            if(MPI.COMM_WORLD.Size() > 1)
            {
                if(rank % 2 == 0)
                {
                    
                    //Even-ranked machines send to the left and right as 
                    //applicable first
                    
                    //if you're not the last machine, send to the right
                    if(rank != MPI.COMM_WORLD.Size() - 1)
                    {
                        //converts last column to one-dimensional 
                        workingColumn = convertColumn(endPos, size, z, p);

                        //sends last column to machine to the right
                        MPI.COMM_WORLD.Send(workingColumn, 0, size, MPI.DOUBLE, 
                                rank + 1, id);
                        //System.out.println("sent to machine " + (rank + 1));
                    }
                    //if you're not first machine, send to the left
                    if(rank != 0)
                    {
                        //converts first column to one-dimensional array
                        workingColumn = convertColumn(startPos, size, z, p);

                        //sends first column to machine to the left
                        MPI.COMM_WORLD.Send(workingColumn, 0, size, MPI.DOUBLE,
                                rank - 1, id);
                        //System.out.println("sent to machine " + (rank - 1));
                    }
                    
                    //Even-ranked machines receive from left and right as
                    //applicable second
                    
                    //if you're not the last machine, receive from the right
                    if(rank != MPI.COMM_WORLD.Size() - 1)
                    {
                        //receives right machine's first column
                        MPI.COMM_WORLD.Recv(workingColumn, 0, size, MPI.DOUBLE, 
                                rank + 1, id);

                        //sets received column in machine's z[][][]
                        setColumn(endPos + 1, size, z, p, workingColumn);

                    }
                    //if you're not the first machine, receive from the left
                    if(rank != 0)
                    {
                        //receives left machine's last column
                        MPI.COMM_WORLD.Recv(workingColumn, 0, size, MPI.DOUBLE, 
                                rank - 1, id);

                        //sets received column in machine's z[][][]
                        setColumn(startPos - 1, size, z, p, workingColumn);
                    }
                }
                
                else
                {
                    //Odd-ranked machines receive from the left and right as 
                    //applicable first
                    
                    //if you're not the last machine, receive from the right
                    if(rank != MPI.COMM_WORLD.Size() - 1)
                    {
                        //receives right machine's first column
                        MPI.COMM_WORLD.Recv(workingColumn, 0, size, MPI.DOUBLE, 
                                rank + 1, id);
                        
                        //sets received column in machine's z[][][]
                        setColumn(endPos + 1, size, z, p, workingColumn);
                    }
                    //if you're not the first machine, receive from the left
                    if(rank != 0)
                    {
                        //receives left machine's last column
                        MPI.COMM_WORLD.Recv(workingColumn, 0, size, MPI.DOUBLE,
                                rank - 1, id);

                        //sets received column in machine's z[][][]
                        setColumn(startPos - 1, size, z, p, workingColumn);

                    }
                    //Odd-ranked machines send to the left and right as 
                    //applicable second
                    
                    //if you're not the last machine, send to the right
                    if(rank != MPI.COMM_WORLD.Size() - 1)
                    {
                        //converts last column to one-dimensional 
                        workingColumn = convertColumn(endPos, size, z, p);

                        //sends last column to machine to the right
                        MPI.COMM_WORLD.Send(workingColumn, 0, size, MPI.DOUBLE,
                                rank + 1, id);
                    }
                    //if you're not first machine, send to the left
                    if(rank != 0)
                    {
                        //converts first column to one-dimensional array
                        workingColumn = convertColumn(startPos, size, z, p);

                        //sends first column to machine to the left
                        MPI.COMM_WORLD.Send(workingColumn, 0, size, MPI.DOUBLE,
                                rank - 1, id);
                    }
                }
            }
        }
        catch(MPIException e)
        {
            System.out.println("An error occured in shareBoundaryData()");
        }
    }




    /*
        setColumn
        takes a one-dimensional array of doubles and changes
        z[p][pos][<variable>] to it
    */
    public static void setColumn(int pos, int size, double[][][] z, int p,
            final double[] column)
    {
        for(int i = 0; i < size; i++)
        {
            z[p][pos][i] = column[i];
        }
    }




    /*
        convertColumn
        Converts an entire column of one of the layers of the three dimensional
        array into a one-dimensional array. This array represents this column
        in a top-down fashion
    
    */
    public static double[] convertColumn(int pos, int size, 
            final double[][][] z, int p)
    {
        double[] workingColumn = new double[size];
        for(int i = 0; i < size; i++)
        {
            workingColumn[i] = z[p][pos][i];
        }
        return workingColumn;
    }




    /*
        determinePositions
        Determines the start and end positions of a square that a machine is
        responsible for
        Used by the master machine, who then communicates this information to
        other computers
    */
    public static void determinePositions(int size)
    {
        //System.out.println("Master entered determinePositions");
        try
        {
            int numOfMachines = MPI.COMM_WORLD.Size();
            int remainder = size % numOfMachines;
            Boolean remainderPresent = (remainder != 0);
            int partition = size / numOfMachines;
            int start = 0;
            int end = 0;
            
            for(int i = 0; i < numOfMachines; i++)
            {
                end += partition - 1;
                
                //if size of square is not evenly divisible into numOfMachines,
                //allocate the remainder to the last machines
                if(i == numOfMachines - 1 && remainderPresent)
                {
                    end = size - 1;
                }
                
                int[] positions = new int[]{start, end};
                
                //sends position information to appropriate machine
                if(i == 0)
                {
                    //do not send to self
                    startPos = positions[0];
                    endPos = positions[1];
                }
                else
                {
                    //send position information to other machines
                    MPI.COMM_WORLD.Send(positions, 0, 2, MPI.INT, i, id);
                }
                
                //sets start and end to be 1 further
                start = end + 1;
                end++;
            } 
        }
        catch(MPIException e)
        {
            System.out.println("An exception occured in determinePosition");
        }
    }
}
