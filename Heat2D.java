import java.util.Date;

public class Heat2D
{
    private static double a = 1.0;  // heat speed
    private static double dt = 1.0; // time quantum
    private static double dd = 2.0; // change in system
    
    public static void main( String[] args )
    {
	// verify arguments
	if ( args.length != 4 )
        {
	    System.out.println( "usage: " + 
			 "java Heat2D size max_time heat_time interval" );
	    System.exit( -1 );
	}

	int size = Integer.parseInt( args[0] );
	int max_time = Integer.parseInt( args[1] );
	int heat_time = Integer.parseInt( args[2] );
	int interval  = Integer.parseInt( args[3] );
	double r = a * dt / ( dd * dd );

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
	
	
	// start a timer
	Date startTime = new Date( );
	
	// simulate heat diffusion
	for ( int t = 0; t < max_time; t++ )
        {
	    int p = t % 2; // p = 0 or 1: indicates the phase
	    
	    // two left-most and two right-most columns are identical
	    for ( int y = 0; y < size; y++ )
            {
		z[p][0][y] = z[p][1][y];
		z[p][size - 1][y] = z[p][size - 2][y];
	    }
	    
	    // two upper and lower rows are identical
	    for ( int x = 0; x < size; x++ )
            {
		z[p][x][0] = z[p][x][1];
		z[p][x][size - 1] = z[p][x][size - 2];
	    }
	    
	    // keep heating the bottom until t < heat_time
	    if ( t < heat_time )
            {
		for ( int x = size /3; x < size / 3 * 2; x++ )
                {
                    z[p][x][0] = 19.0; // heat
                }
            }

	    // display intermediate results
	    if ( interval != 0 && 
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
	    
	    // perform forward Euler method
	    int p2 = (p + 1) % 2;
	    for ( int x = 1; x < size - 1; x++ )
            {
                for ( int y = 1; y < size - 1; y++ )
                {
                    z[p2][x][y] = z[p][x][y] + 
			r * ( z[p][x+1][y] - 2 * z[p][x][y] + z[p][x-1][y] ) +
			r * ( z[p][x][y+1] - 2 * z[p][x][y] + z[p][x][y-1] );
                }  
            }
	 
	} // end of simulation
	
	// finish the timer
	Date endTime = new Date( );
	System.out.println( "Elapsed time = " + 
			    ( endTime.getTime( ) - startTime.getTime( ) ) );
    }
}
