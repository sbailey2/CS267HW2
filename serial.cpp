#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//

extern double size;

#define NUM_BLOCKS 20

int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    particle_t **first = (particle_t**) malloc( NUM_BLOCKS * NUM_BLOCKS * sizeof(particle_t*) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
	organize_particles(n, particles, first, NUM_BLOCKS);
        //
        //  compute forces
        //
        for( int i = 0; i < NUM_BLOCKS * NUM_BLOCKS; i++ )
        {
	    unsigned int x = i / NUM_BLOCKS;
	    unsigned int y = i % NUM_BLOCKS;
	    for (particle_t *cur = first[i]; cur != 0; cur = cur->next) {
		cur->ax = cur->ay = 0;

		// Apply forces on particles in the grid
		for (particle_t *neighbor = first[i]; neighbor != 0; neighbor = neighbor->next) {
		    apply_force(*cur, *neighbor, &dmin, &davg, &navg);
		}
	    }
	  
	    // Apply forces on particles on the border
	    particle_t **neighbors[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	    unsigned int n[8] = {0, 0, 0, 0, 0, 0, 0};
	    // Upper left
	    if (x > 0 && y > 0) {
		neighbors[0] = get_neighboring_particles(first[x - 1 + (y - 1) * NUM_BLOCKS], DOWN | RIGHT, n[0]);
	    }
	    // Up
	    if (y > 0) {
		neighbors[1] = get_neighboring_particles(first[x + (y - 1) * NUM_BLOCKS], DOWN, n[1]);
	    }
	    // Upper right
	    if (x < NUM_BLOCKS - 1 && y > 0) {
		neighbors[2] = get_neighboring_particles(first[x + 1 + (y - 1) * NUM_BLOCKS], DOWN | LEFT, n[2]);
	    }
	    // LEFT
	    if (x < 0) {
		neighbors[3] = get_neighboring_particles(first[x - 1 + y * NUM_BLOCKS], RIGHT, n[3]);
	    }
	    // RIGHT
	    if (x < NUM_BLOCKS - 1) {
		neighbors[4] = get_neighboring_particles(first[x + 1 + y * NUM_BLOCKS], LEFT, n[4]);
	    }
	    // Lower left
	    if (x > 0 && y < NUM_BLOCKS - 1) {
		neighbors[5] = get_neighboring_particles(first[x - 1 + (y + 1) * NUM_BLOCKS], UP | RIGHT, n[5]);
	    }
	    // Down
	    if (y < NUM_BLOCKS - 1) {
		neighbors[6] = get_neighboring_particles(first[x + (y + 1) * NUM_BLOCKS], UP, n[6]);
	    }
	    // Lower right
	    if (x < NUM_BLOCKS - 1 && y < NUM_BLOCKS - 1) {
		neighbors[7] = get_neighboring_particles(first[x + 1 + (y + 1) * NUM_BLOCKS], UP | LEFT, n[7]);
	    }

	    // Calculate forces from neighboring grid cells
	    for (particle_t *cur = first[i]; cur != 0; cur = cur->next) {
		if (cur->boundary & LEFT != 0) {
		    for (unsigned int j = 0; j < n[3]; ++j) {
			apply_force(*cur, *neighbors[3][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & RIGHT != 0) {
		    for (unsigned int j = 0; j < n[4]; ++j) {
			apply_force(*cur, *neighbors[4][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & UP != 0) {
		    for (unsigned int j = 0; j < n[1]; ++j) {
			apply_force(*cur, *neighbors[1][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & DOWN != 0) {
		    for (unsigned int j = 0; j < n[6]; ++j) {
			apply_force(*cur, *neighbors[6][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & UP != 0 && cur->boundary & LEFT != 0) {
		    for (unsigned int j = 0; j < n[0]; ++j) {
			apply_force(*cur, *neighbors[0][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & UP != 0 && cur->boundary & RIGHT != 0) {
		    for (unsigned int j = 0; j < n[2]; ++j) {
			apply_force(*cur, *neighbors[2][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & DOWN != 0 && cur->boundary & LEFT != 0) {
		    for (unsigned int j = 0; j < n[5]; ++j) {
			apply_force(*cur, *neighbors[5][j], &dmin, &davg, &navg);
		    }
		}
		if (cur->boundary & DOWN != 0, cur->boundary & RIGHT != 0) {
		    for (unsigned int j = 0; j < n[7]; ++j) {
			apply_force(*cur, *neighbors[7][j], &dmin, &davg, &navg);
		    }
		}
	    }

	    // Clean up the neighbors
	    for (unsigned int i = 0; i < 8; ++i) {
		free(neighbors[i]);
	    }
	  
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( first );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
