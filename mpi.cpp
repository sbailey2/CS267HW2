#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include "common.h"

using std::vector;

extern double size;
MPI_Datatype PARTICLE;

int M, N;
double pRowWidth, pColWidth;
double nGBdry, sGBdry, eGBdry, wGBdry;

void find_dim(double size, int n_proc) {
    M = N = 1;
    bool multN = 1;
    while (n_proc != 1) {
	for (int i = 2; i <= n_proc; i++) {
	    if (n_proc % i == 0) {
		n_proc = n_proc / i;
		(multN) ? (N=N*i) : (M=M*i);
		multN = 1-multN;	
		break;
	    }
	}
    }
}

void find_neighbor_procs(int i, int j, int *pNeighbors) {
    pNeighbors[0] = ((i>0)&&(j>0)) ? ((i-1)*N+(j-1)) : (-1);
    pNeighbors[1] = (i>0) ? ((i-1)*N+j) : (-1);
    pNeighbors[2] = ((i>0)&&(j<N-1)) ? ((i-1)*N+(j+1)) : (-1);
    pNeighbors[3] = (j>0) ? (i*N+(j-1)) : (-1);
    pNeighbors[4] = (j<N-1) ? (i*N+(j+1)) : (-1);
    pNeighbors[5] = ((i<M-1)&&(j>0)) ? ((i+1)*N+(j-1)) : (-1);
    pNeighbors[6] = (i<M-1) ? ((i+1)*N+j) : (-1);
    pNeighbors[7] = ((i<M-1)&&(j<N-1)) ? ((i+1)*N+(j+1)) : (-1);
}

void setup_proc_grid(double size, int n_proc, int rank, int *pRowInd, int *pColInd, int *pNeighbors) {

    find_dim(size, n_proc);
    (*pRowInd) = rank/N;
    (*pColInd) = rank - (*pRowInd)*N;
    pRowWidth = size/(double)M;
    pColWidth = size/(double)N;
    find_neighbor_procs(*pRowInd,*pColInd,pNeighbors);
    nGBdry = ((*pRowInd)+1)*(pRowWidth) - cutoff;
    sGBdry = (*pRowInd)*(pRowWidth) + cutoff;
    eGBdry = ((*pColInd)+1)*(pColWidth) - cutoff;
    wGBdry = (*pColInd)*(pColWidth) + cutoff;
}

int get_proc_index(particle_t &p) {
    int i = p.y/pRowWidth;
    int j = p.x/pColWidth;
    return i*N+j;
}

void init_map_particles_to_procs(vector<vector<particle_t> > *dist, particle_t *particles, int n) {
 
    dist->resize(M*N);
    for (int i = 0; i < n; i++) {
	int j = get_proc_index(particles[i]);
	(*dist)[j].push_back(particles[i]);
    }
}

// debugging function
void print_2DGrid(vector<vector<particle_t> > dist) {

    int total = 0;
    for (int i = 0; i < dist.size(); i++) {
	printf("proc %d will have %d particles.\n",i,dist[i].size());
        total += dist[i].size();
    }
    printf("total particles: %d\n",total);
}

void find_ghostzone(vector<vector<particle_t> > *ghostzone, vector<particle_t> *local, int n) {

    int x,y;
    for (int i = 0; i < (*local).size(); i++) {
	x = (*local)[i].x;
	y = (*local)[i].y;
	if ((x < wGBdry) && (y < sGBdry))
	    (*ghostzone)[0].push_back((*local)[i]);
	if (y < sGBdry)
	    (*ghostzone)[1].push_back((*local)[i]);
	if ((x > eGBdry) && (y < sGBdry))
	    (*ghostzone)[2].push_back((*local)[i]);
	if (x < wGBdry)
	    (*ghostzone)[3].push_back((*local)[i]);
	if (x > eGBdry)
	    (*ghostzone)[4].push_back((*local)[i]);
	if ((x < wGBdry) && (y > nGBdry))
	    (*ghostzone)[5].push_back((*local)[i]);
	if (y > nGBdry)
	    (*ghostzone)[6].push_back((*local)[i]);
	if ((x > eGBdry) && (y > nGBdry))
	    (*ghostzone)[7].push_back((*local)[i]);
    }
}

int pick_ind(int r, int c, int i) {
    if ((i==1) || (i==6)) {
	return r;
    } else {
	return c;
    }
}

// THIS FUNCTION IS NOT WORKING AND NEEDS TO BE DEBUGGED
void send_and_recv_ghostzone(vector<vector<particle_t> > *out, vector<vector<particle_t> > *in, int* neighbors, int pRowInd, int pColInd,int n) {
    
    int nSendSize;
    int di;
    int ind;
    for (int i = 1; i < 2; i++) {
    // even columns send particles, old columns receive particles
    di = 7-i;
    ind = pick_ind(pRowInd, pColInd,i);
    if (((ind%2)==0) && (neighbors[i]!=-1)) {
	printf("proc (%d,%d) sending to %d \n",pRowInd,pColInd,neighbors[i]);
	nSendSize = (*out)[i].size();
	MPI_Send(&nSendSize,1,MPI_UNSIGNED,neighbors[i],0,MPI_COMM_WORLD);
	MPI_Send(&out[i].front(),nSendSize,PARTICLE,neighbors[i],0,MPI_COMM_WORLD);
    }
    if (((ind%2)!=0) && (neighbors[di]!=-1)) {
	printf("proc (%d,%d) receiving from %d\n",pRowInd,pColInd,neighbors[di]);
	MPI_Recv(&nSendSize,1,MPI_UNSIGNED,neighbors[di],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	(*in)[di].resize(nSendSize);
	MPI_Recv(&in[di].front(),nSendSize,PARTICLE,neighbors[di],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    // odd columns send particles, even columns receive particles
/*
    if (((pColInd%2)!=0) && (neighbors[i]!=-1)) {
	nSendSize = (*out)[i].size();
	MPI_Send(&nSendSize,1,MPI_UNSIGNED,neighbors[i],0,MPI_COMM_WORLD);
	MPI_Send(&out[i],nSendSize,PARTICLE,neighbors[i],0,MPI_COMM_WORLD);
    }
    if (((pColInd%2)==0) && (neighbors[dual(i)]!=-1)) {
	MPI_Recv(&nSendSize,1,MPI_UNSIGNED,neighbors[dual(i)],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	(*in)[dual(i)].resize(nSendSize);
	MPI_Recv(&in[dual(i)].front(),n,PARTICLE,neighbors[dual(i)],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    } */
    }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;
     
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    particle_t *particles;
    vector<vector<particle_t> > particle_dist;
    vector<particle_t> local;
    vector<vector<particle_t> > ghostzone(8);
    vector<vector<particle_t> > incoming_ghostzone(8);
    
    // 
    // set up partitioning across processors based on grid size
    //
    int pRowInd, pColInd;
    int pNeighbors[8];
    set_size(n);
    setup_proc_grid(size, n_proc, rank, &pRowInd, &pColInd, pNeighbors);
    
    // debugging sanity check
    //printf("Hi, I'm processor %d of %d. In the %d x %d processing grid, my coordinates are (%d,%d). The grid is %f x %f, and each proc is in charge of an area of %f x %f. My neighbors are: %d,%d,%d,%d,%d,%d,%d,%d.\n ",rank,n_proc,M,N,pRowInd,pColInd,size,size,pRowWidth,pColWidth,pNeighbors[0],pNeighbors[1],pNeighbors[2],pNeighbors[3],pNeighbors[4],pNeighbors[5],pNeighbors[6],pNeighbors[7]);

    // initialize particles on proc 0 and send particles to each of the other procs
    if (rank == 0) {
	particles = (particle_t*) malloc( n * sizeof(particle_t) );
        init_particles(n, particles);
        init_map_particles_to_procs(&particle_dist, particles, n);
	local = particle_dist[0];
//	print_2DGrid(particle_dist);

	for (int i = 1; i < n_proc; i++) {
	    int nSendSize = particle_dist[i].size();
	    MPI_Send(&nSendSize,1,MPI_UNSIGNED,i,0,MPI_COMM_WORLD); // send size
	    MPI_Send(&particle_dist[i].front(),particle_dist[i].size(),PARTICLE,i,0,MPI_COMM_WORLD); // send particles
	}
    } else {
	int nlocal;
	MPI_Recv(&nlocal,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // receive size
	local.resize(nlocal); // resize to right length
	MPI_Recv(&local.front(),n,PARTICLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // receive particles
    }    

//    printf("Proc %d has %d particles.\n",rank,local.size());

    //
    //  set up the data partitioning across processors
    //
    /*
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    */

    //
    //  allocate storage for local partition
    //
    /*
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    */

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    /*
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    */

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < 1; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        // 
        //  collect all global data locally (not good idea to do)
        //
        //MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
//        if( find_option( argc, argv, "-no" ) == -1 )
//          if( fsave && (step%SAVEFREQ) == 0 )
//            save( fsave, n, particles );


	// 
	//  find ghostzone particles and send to the right processors
	//

	find_ghostzone(&ghostzone,&local,n);	



//	send_and_recv_ghostzone(&ghostzone,&incoming_ghostzone,pNeighbors,pRowInd,pColInd,n);
	
        
        //
        //  compute all forces
        //
//        for( int i = 0; i < nlocal; i++ )
//        {
//            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++ )
//                apply_force( local[i], particles[j], &dmin, &davg, &navg );
//        }
     
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        //
        //  move particles
        //
//        for( int i = 0; i < nlocal; i++ )
//            move( local[i] );
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    //free( partition_offsets );
    //free( partition_sizes );
    //free( local );
    if (rank == 0)
	free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
