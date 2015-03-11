#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include "common2.h"

using std::vector;

extern double size;
int M, N;
int pRowInd, pColInd;
double pRowWidth, pColWidth;
double nBdry, sBdry, eBdry, wBdry;
double ngBdry, sgBdry, egBdry, wgBdry;
int pNeighbors[8];
int nn;

int gridm, gridn;
double gridColWidth, gridRowWidth;
int offset[9];

void setup_local_grid() {
    gridm = (int)(pRowWidth/cutoff);
    gridn = (int)(pColWidth/cutoff);
    gridColWidth = pColWidth/(double)gridn;
    gridRowWidth = pRowWidth/(double)gridm;
    gridm += 2;
    gridn += 2;
}

int grid_index(double x, double y) {
    int cInd = (int)((x-wBdry)/gridColWidth)+1;
    int rInd = (int)((y-sBdry)/gridRowWidth)+1;
    return rInd*gridn+cInd;
}

void compute_offset() {

    offset[0] = (-1)*gridn+(-1);
    offset[1] = (-1)*gridn;
    offset[2] = (-1)*gridn+(1);
    offset[3] = -1;
    offset[4] = 0;
    offset[5] = 1;
    offset[6] = (1)*gridn+(-1);
    offset[7] = (1)*gridn;
    offset[8] = (1)*gridn+(1);
}

void local_apply_force(particle_t &self, vector<particle_t> &neighbor, double *dmin, double *davg, int *navg) {
    for (int j = 0; j < neighbor.size(); j++) {
	apply_force(self,neighbor[j],dmin,davg,navg);
    }
}

void distribute_particles_to_local_grid(vector<particle_t> local, vector<vector<particle_t> > gz, vector<vector<particle_t> > &grid) {

    // clear grid from last time
    for (int i = 0; i < gridm*gridn; i++) {
	grid[i].clear();
    }

    // bin particles from local and ghostzone into grid
    for (int i = 0; i < local.size(); i++) {
	grid[grid_index(local[i].x,local[i].y)].push_back(local[i]);
    }
    vector<particle_t> currGZ;
    for (int i = 0; i < 8; i++) {
	currGZ = gz[i];
	for (int j = 0; j < currGZ.size(); j++) {
	    grid[grid_index(currGZ[j].x,currGZ[j].y)].push_back(currGZ[j]);
	}
    }
}

void find_dim(int n_proc) {
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

void find_neighbor_procs() {
    pNeighbors[0] = ((pRowInd>0)&&(pColInd>0)) ? ((pRowInd-1)*N+(pColInd-1)) : (-1);
    pNeighbors[1] = (pRowInd>0) ? ((pRowInd-1)*N+pColInd) : (-1);
    pNeighbors[2] = ((pRowInd>0)&&(pColInd<N-1)) ? ((pRowInd-1)*N+(pColInd+1)) : (-1);
    pNeighbors[3] = (pColInd>0) ? (pRowInd*N+(pColInd-1)) : (-1);
    pNeighbors[4] = (pColInd<N-1) ? (pRowInd*N+(pColInd+1)) : (-1);
    pNeighbors[5] = ((pRowInd<M-1)&&(pColInd>0)) ? ((pRowInd+1)*N+(pColInd-1)) : (-1);
    pNeighbors[6] = (pRowInd<M-1) ? ((pRowInd+1)*N+pColInd) : (-1);
    pNeighbors[7] = ((pRowInd<M-1)&&(pColInd<N-1)) ? ((pRowInd+1)*N+(pColInd+1)) : (-1);
    nn = 0;
    for (int i = 0; i < 8; i++) 
	nn += (pNeighbors[i]!=-1) ? (1) : (0);
}

void setup(int n_proc, int rank) {

    find_dim(n_proc);
    pRowInd = rank/N;
    pColInd = rank - pRowInd*N;
    pRowWidth = size/(double)M;
    pColWidth = size/(double)N;
    find_neighbor_procs();
    nBdry = (pRowInd+1)*pRowWidth;
    sBdry = pRowInd*pRowWidth;
    eBdry = (pColInd+1)*pColWidth;
    wBdry = pColInd*pColWidth;
    ngBdry = nBdry - cutoff;
    sgBdry = sBdry + cutoff;
    egBdry = eBdry - cutoff;
    wgBdry = wBdry + cutoff;
}

void map_particles_to_processors(int n_proc, int n, particle_t *in, int *sizes, int *offsets, vector<particle_t> &out) {
    vector<vector<particle_t> > dist(n_proc);
    int rInd, cInd;
    for (int i = 0; i < n; i++) {
	cInd = in[i].x/pColWidth;
	rInd = in[i].y/pRowWidth;
	dist[rInd*N+cInd].push_back(in[i]);
    }
    offsets[0]=0;
    for (int i = 0; i < n_proc; i++) {
	vector<particle_t> curr = dist[i];
	for (int j = 0; j < curr.size(); j++) {
	    out.push_back(curr[j]);
	}
	if (i != n_proc-1)
	    offsets[i+1] = offsets[i]+curr.size();
	sizes[i] = curr.size();
    }
}

void compute_particle_distribution(vector<particle_t> &send, int *dist_sizes, int *dist_disp, vector<particle_t> local, int n_proc) {

    int rInd, cInd;
    vector<vector<particle_t> > dist(n_proc);
    vector<particle_t> currList;
    send.clear();
    for (int i = 0; i < local.size(); i++) {
	cInd = local[i].x/pColWidth;
	rInd = local[i].y/pRowWidth;
	dist[rInd*N+cInd].push_back(local[i]);
    }
    dist_disp[0]=0;
    for (int i = 0; i < n_proc; i++) {
	currList = dist[i];
	for (int j = 0; j < currList.size(); j++) {
	    send.push_back(currList[j]);
	}
	dist_sizes[i] = currList.size();
	if (i != n_proc-1)
	    dist_disp[i+1] = dist_disp[i]+currList.size();
    }
} 

void compute_particle_gather(vector<particle_t> &recv, int *dist_sizes, int *dist_disp, int n_proc) {
    
    dist_disp[0]=0;
    for (int i = 1; i < n_proc; i++) {
	dist_disp[i] = dist_disp[i-1]+dist_sizes[i-1];
    }
    recv.clear();
    recv.resize(dist_disp[n_proc-1]+dist_sizes[n_proc-1]);
}

void find_ghostzone(vector<vector<particle_t> > &ghostzone, vector<particle_t> local) {

    for (int i = 0; i < 8; i++) {
	ghostzone[i].clear();	
    }
    for (int i = 0; i < local.size(); i++) {
	if (local[i].x > egBdry)
	    ghostzone[4].push_back(local[i]);
	if (local[i].x < wgBdry)
	    ghostzone[3].push_back(local[i]);
	if (local[i].y > ngBdry)
	    ghostzone[6].push_back(local[i]);
	if (local[i].y < sgBdry)
	    ghostzone[1].push_back(local[i]);
	if ((local[i].x < wgBdry) && (local[i].y < sgBdry))
	    ghostzone[0].push_back(local[i]);
	if ((local[i].x > egBdry) && (local[i].y < sgBdry))
	    ghostzone[2].push_back(local[i]);
	if ((local[i].x < wgBdry) && (local[i].y > ngBdry))
	    ghostzone[5].push_back(local[i]);
	if ((local[i].x > egBdry) && (local[i].y > ngBdry))
	    ghostzone[7].push_back(local[i]);
    }
}

void generate_disp(int *count, int n, int *disp) {
    disp[0] = 0;
    for (int i = 1; i < n; i++) {
	disp[i] = disp[i-1]+count[i-1];
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


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
   

    //  
    //  set up constants and rank dependent parameters
    //
    set_size(n);
    setup(n_proc, rank);
    setup_local_grid();
    compute_offset();

    //
    //  initialize variables 
    //
    int partition_offsets[n_proc];
    int partition_sizes[n_proc];
    vector<particle_t> particles_reordered; 
    vector<particle_t> local;
    vector<particle_t> dist_send;
    vector<particle_t> dist_recv;
    vector<particle_t> currNeighbor;
    vector<vector<particle_t> > iGZ(8);
    vector<vector<particle_t> > oGZ(8);
    vector<vector<particle_t> > local_grid(gridm*gridn);
    int dist_send_count[n_proc];
    int dist_recv_count[n_proc];
    int dist_send_disp[n_proc];
    int dist_recv_disp[n_proc];
    int nlocal;
    int save_count[n_proc];
    int save_disp[n_proc];
    int saveN;

    //
    //  distribute particles to each processor
    //
    if (rank == 0) {
	init_particles(n,particles);
	map_particles_to_processors(n_proc,n,particles,partition_sizes,partition_offsets,particles_reordered);
    }
    MPI_Scatter(&partition_sizes[0],1,MPI_INT,&nlocal,1,MPI_INT,0,MPI_COMM_WORLD);
    local.resize(nlocal);
    MPI_Scatterv(&particles_reordered.front(), &partition_sizes[0], &partition_offsets[0], PARTICLE, &local.front(), nlocal, PARTICLE, 0, MPI_COMM_WORLD);
  
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
       
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 ) {
          if( savename && (step%SAVEFREQ) == 0 ) {
	    saveN = local.size();
	    MPI_Gather(&saveN,1,MPI_INT,&save_count[0],1,MPI_INT,0,MPI_COMM_WORLD);
	    if (rank == 0) { 
    	        generate_disp(save_count,n_proc,save_disp);
	    }
	    MPI_Gatherv(&local.front(),local.size(),PARTICLE,particles,save_count,save_disp,PARTICLE,0,MPI_COMM_WORLD);
  	    if (rank == 0) {
		save( fsave, n, particles );
	    }
	  }
	}

	// 
	//  calculate ghostzone for each processor and distribute
	//
	find_ghostzone(oGZ,local);

	MPI_Request reqs[2*nn];
	MPI_Status stats[2*nn];
	int counter = 0;
	for (int i = 0; i < 8; i++) {
	    if (pNeighbors[i]!=-1) {
		MPI_Isend(&oGZ[i].front(),oGZ[i].size(),PARTICLE,pNeighbors[i],0,MPI_COMM_WORLD,&reqs[counter]);
		counter++;
		iGZ[i].resize(n);
		MPI_Irecv(&iGZ[i].front(),n,PARTICLE,pNeighbors[i],0,MPI_COMM_WORLD,&reqs[counter]);
		counter++;
	    }
	}
	MPI_Waitall(2*nn,reqs,stats);        
	int np;
	counter = 1;
	for (int i = 0; i < 8; i++) {
	    if (pNeighbors[i]!=-1) {
		MPI_Get_count(&stats[counter],PARTICLE,&np);
		iGZ[i].resize(np);
		counter+=2;
	    }
	}	

        //
        //  compute all forces using nlocal interactions
        //
	
	int gInd;
	distribute_particles_to_local_grid(local,iGZ,local_grid);
	for (int i = 0; i < local.size(); i++) {
	    local[i].ax = local[i].ay = 0;
	    gInd = grid_index(local[i].x,local[i].y);
	    local_apply_force(local[i],local_grid[gInd+offset[0]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[1]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[2]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[3]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[4]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[5]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[6]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[7]],&dmin,&davg,&navg);
	    local_apply_force(local[i],local_grid[gInd+offset[8]],&dmin,&davg,&navg);
	}

	//
	//  compute all forces naively using nlocal^2 interactions
	//
/*	 
        for( int i = 0; i < local.size(); i++ )
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < local.size(); j++ )
                apply_force( local[i], local[j], &dmin, &davg, &navg );
	    for (int k = 0; k < 8; k++) {
		currNeighbor = iGZ[k];
		for (int j = 0; j < currNeighbor.size(); j++) 
		    apply_force(local[i],currNeighbor[j],&dmin,&davg,&navg);
	    }
        }
	
*/
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
        for( int i = 0; i < local.size(); i++ )
            move( local[i] );

	//
	//  send particles to the right processor after they move
	//
	compute_particle_distribution(dist_send,dist_send_count,dist_send_disp,local,n_proc);
	MPI_Alltoall(&dist_send_count[0],1,MPI_INT,&dist_recv_count[0],1,MPI_INT,MPI_COMM_WORLD);
	compute_particle_gather(dist_recv,dist_recv_count,dist_recv_disp,n_proc);
	MPI_Alltoallv(&dist_send.front(),dist_send_count,dist_send_disp,PARTICLE,&dist_recv.front(),dist_recv_count,dist_recv_disp,PARTICLE,MPI_COMM_WORLD);
	local = dist_recv;
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
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
