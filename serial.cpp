#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>

// function to place particles into the correct blocks
// n = number of particles
// N = number of blocks in one row or column
// radius = length and width of block
void bin_particles(int n, int N, double radius, particle_t *particles, std::vector<std::vector<int>> &blocks)
{
    for (std::vector<std::vector<int>>::iterator it = blocks.begin(); it != blocks.end(); it++)
    {
	(*it).clear();
    }
    int rowind, colind;
    particle_t particle;
    for (int i = 0; i < n; ++i) {
        particle = particles[i];
        rowind = (int) (particle.x/radius);
        colind = (int) (particle.y/radius);
	blocks[colind+rowind*N].push_back(i);
    }
}

// debugging helper function
// verifies that the right number of particles are present
void verify_bins(std::vector<std::vector<int>> &blocks, int n)
{
    int total = 0;
    for (std::vector<std::vector<int>>::iterator it = blocks.begin(); it != blocks.end(); it++)
    {
	total += (*it).size();
    }
    if (total != n)
	printf("You have the wrong number of particles. Have: %d; Want: %d.", total, n);
}

// debugging helper function
void print_bins(std::vector<std::vector<int>> &blocks,int N)
{
    int particle;
    int total = 0;
    int counter = 0;
    std::vector<int> blockrow;
    for (std::vector<std::vector<int>>::iterator it = blocks.begin(); it != blocks.end(); ++it) {
	int counter2 = 0;        
	blockrow = *it;
	for (std::vector<int>::iterator it2 = blockrow.begin(); it2 != blockrow.end(); ++it2) {
	    particle = (int) *it2;	    
	    printf("block %d: %d\n",counter,particle);
	    counter2++;
	}
	printf("block %d has %d particles\n",counter,counter2);
	total+=counter2;
	counter++;
    } 
    printf("there is a total of %d blocks\n",counter);
    printf("there is a total of %d particles\n",total);
} 

// function to let particles in self and neighbor interact
void interact_particles(particle_t *particles, std::vector<int> &self, std::vector<int> &neighbor, double &dmin, double &davg, int &navg)
{
    int pi, pj;
    
    for (std::vector<int>::iterator it1 = self.begin(); it1 != self.end(); it1++)
    {
	pi =  *it1;
	for (std::vector<int>::iterator it2 = neighbor.begin(); it2 !=neighbor.end(); it2++)
	{
	    pj =  *it2;
	    //printf("(%d,%d)",pi,pj);
	    apply_force(particles[pi],particles[pj],&dmin,&davg,&navg);
	}
    }
}

void check_simulation(particle_t *particles,int n,double &dmin,double &davg,int &navg)
{

        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }

}

//
//  benchmarking program
//
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
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    const double size = get_size();
    const int N = (int) (size/cutoff); // grid is split into N x N blocks
   const double radius = size/(double)N;
//	printf("radius %f, cutoff %f\n",radius,cutoff);

//const int N = 1;
//const double radius = size;
    // initialize N*N blocks    
    std::vector<std::vector<int>> blocks(N*N);
    bin_particles(n,N,radius,particles,blocks);
    verify_bins(blocks, n);

    // debugging, check that particles are put into bins correctly
//     print_bins(blocks,N);

    std::vector<int> self, no, so, e, w, ne, nw, se, sw;
    int i,j;

    //printf("%d\n",N);

    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces by looping through blocks
        //

	// first, you need to reset the forces on all the particles
	for (int i = 0; i < n; i++ )
	{
	    particles[i].ax = particles[i].ay = 0;
	}

	for (std::vector<std::vector<int>>::size_type it = 0; it != (N*N); it++)
	{
	    self = blocks[it];
	    if (self.size() > 0)
	    {
	        i = (int) (it/N);
	        j = it - i*N;
	  	//printf("iteration %d: (%d, %d)\n",(int)it,i,j);
		
		// interact with self
		interact_particles(particles,self,self,dmin,davg,navg);
	    
		// interact with northern neighbors
	    	if (i < N-1) 
	    	{
		    no = blocks[j+(i+1)*N];
		    if (no.size() > 0) interact_particles(particles,self,no,dmin,davg,navg);
		    
		    // interact with eastern neighbors
		    if (j < N-1)
		    {
		    	ne = blocks[j+1+(i+1)*N];
		    	if (ne.size() > 0) interact_particles(particles,self,ne,dmin,davg,navg);
		    	e = blocks[j+1+i*N];
		    	if (e.size() > 0) interact_particles(particles,self,e,dmin,davg,navg);
		    }
		    // interact with western neighbors
		    if (j > 1)
		    {
		    	nw = blocks[j-1+(i+1)*N];
		    	if (nw.size() > 0) interact_particles(particles,self,nw,dmin,davg,navg);
		    	w = blocks[j-1+i*N];
		    	if (w.size() > 0) interact_particles(particles,self,w,dmin,davg,navg);
		    }
	    	}
	    	// interact with southern neighbors
	    	if (i > 1)
	    	{
		    so = blocks[j+(i-1)*N];
		    if (so.size() > 0) interact_particles(particles,self,so,dmin,davg,navg);
		
		    // interact with eastern neighbors
	 	    if (j < N-1)
		    {
		    	se = blocks[j+1+(i-1)*N];
		    	if (se.size() > 0) interact_particles(particles,self,se,dmin,davg,navg);
		    }
		    // interact with western neighbors
		    if (j > 1)
		    {
		    	sw = blocks[j-1+(i-1)*N];
		    	if (sw.size() > 0) interact_particles(particles,self,sw,dmin,davg,navg);
		    }
	    	}
	    }
	}

//        for( int i = 0; i < n; i++ )
//        {
//            particles[i].ax = particles[i].ay = 0;
//            for (int j = 0; j < n; j++ )
//				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
//        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

	// rebin after every step
	bin_particles(n,N,radius,particles,blocks);
	verify_bins(blocks,n);

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


// navg,nabsavg=0;
// davg,dmin, absmin=1.0, absavg=0.0;

printf("\n");
//    check_simulation(particles,n,dmin,davg,navg);

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
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
