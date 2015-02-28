#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

double get_size()
{
    return size;
}

// 
// Initialize bins
//
void initialize_bins(bin_t *bins, int N)
{
    int rowind, colind;
    int N2 = N*N;

    for (int k = 0; k < N2; k++)
    {
	rowind = k/N;
	colind = k - rowind*N;
	bins[k].i = rowind;
	bins[k].j = colind;
	bins[k].particles = NULL;
	if (rowind > 0) {
	    bins[k].south = &bins[(rowind-1)*N+colind];
	    if (colind > 0) {
		bins[k].west = &bins[rowind*N+(colind-1)];
		bins[k].sw = &bins[(rowind-1)*N+(colind-1)];
	    } else {
		bins[k].west = NULL;
		bins[k].nw = NULL;
		bins[k].ne = NULL;
	    }
	    if (colind < N-1) {
		bins[k].east = &bins[rowind*N+(colind+1)];
		bins[k].se = &bins[(rowind-1)*N+(colind+1)];
	    } else {
		bins[k].east = NULL;
		bins[k].se = NULL;
		bins[k].ne = NULL;
	    }
	} else {
	    bins[k].south = NULL;
	    bins[k].se = NULL;
	    bins[k].sw = NULL;
	}
	if (rowind < N-1) {
	    bins[k].north = &bins[(rowind+1)*N+colind];
	    if (colind < N-1) {
		bins[k].east = &bins[(rowind*N+(colind+1))];
		bins[k].ne = &bins[(rowind+1)*N+(colind+1)];
	    } else {
		bins[k].east = NULL;
		bins[k].se = NULL;
		bins[k].ne = NULL;
	    }
	    if (colind > 0) {
		bins[k].west = &bins[(rowind*N+(colind-1))];
		bins[k].nw = &bins[(rowind+1)*N+(colind-1)];
	    } else {
		bins[k].west = NULL;
		bins[k].nw = NULL;
		bins[k].sw = NULL;
	    }
	} else {
	    bins[k].north = NULL;
	    bins[k].ne = NULL;
	    bins[k].nw = NULL;
	}
    }
}

// 
//  Put particles into bins depending on their position
//
void bin_particles(bin_t *bins, particle_t *particles, int N, int n, double bin_size)
{
    particle_t *temp;
    int rowind, colind;
    for (int i = 0; i < n; i++)
    {
	colind = (int)(particles[i].x/bin_size);
	rowind = (int)(particles[i].y/bin_size);
	temp = bins[rowind*N+colind].particles;
	bins[rowind*N+colind].particles = &particles[i];
	particles[i].next = temp;
    } 
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
 
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
// interact particles from bin self with particles from bin neighbor
// 
void apply_force_bin(particle_t *self, particle_t *neighbor, bool clearForces, double *dmin, double *davg, int *navg)
{
    if ((self == NULL) || (neighbor == NULL))
	return;

    particle_t *temp = self;
    if (clearForces)
    {
	while (temp != NULL) {
	    (*temp).ax = (*temp).ay = 0;
	    temp = (*temp).next;
	}
    }

    particle_t *selfptr = self;
    particle_t *neighborptr = neighbor;
    while (selfptr != NULL) {
	while (neighborptr != NULL) {
	    apply_force(*selfptr,*neighborptr,dmin,davg,navg);
	    neighborptr = (*neighborptr).next;
	    if (selfptr != neighborptr) {
	    }
	}
	selfptr = (*selfptr).next;
	neighborptr = neighbor;
    }
}

//
//  interact binned particles with neighboring binned particles
//
void interact_binned_particles(bin_t *bins, int N, double *dmin, double *davg, int *navg)
{
    int N2 = N*N;
    for (int i = 0; i < N2; i++)
    {
	apply_force_bin(bins[i].particles,bins[i].particles,1,dmin,davg,navg);
	if (bins[i].north != NULL) 
	    apply_force_bin(bins[i].particles,(*bins[i].north).particles,0,dmin,davg,navg);
	if (bins[i].south != NULL)
	    apply_force_bin(bins[i].particles,(*bins[i].south).particles,0,dmin,davg,navg);
 	if (bins[i].east != NULL) 
	    apply_force_bin(bins[i].particles,(*bins[i].east).particles,0,dmin,davg,navg);
	if (bins[i].west != NULL)
	    apply_force_bin(bins[i].particles,(*bins[i].west).particles,0,dmin,davg,navg);   
 	if (bins[i].ne != NULL) 
	    apply_force_bin(bins[i].particles,(*bins[i].ne).particles,0,dmin,davg,navg);
	if (bins[i].nw != NULL)
	    apply_force_bin(bins[i].particles,(*bins[i].nw).particles,0,dmin,davg,navg);
 	if (bins[i].se != NULL) 
	    apply_force_bin(bins[i].particles,(*bins[i].se).particles,0,dmin,davg,navg);
	if (bins[i].sw != NULL)
	    apply_force_bin(bins[i].particles,(*bins[i].sw).particles,0,dmin,davg,navg);   
    }
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

// 
// check if particle is in the correct bin
//
bool isCorrectlyBinned(double x,double y,int i,int j,double bin_size)
{
    if ((y>(i*bin_size)) && (y<=((i+1)*bin_size)) && (x>(j*bin_size)) && (x<=((j+1)*bin_size))) {
	return 1;
    } else {
	return 0;
    }
}

//
// move each particle, then reassign to the correct bin
//
void move_and_rebin(bin_t *bins, int N, double bin_size)
{
    particle_t *rebinQueue = NULL;
    particle_t *stayQueue = NULL;
    particle_t *temp;
    particle_t *particle;
    int N2 = N*N;
    int i, j;

    // for each bin
    for (int i = 0; i < N2; i++) {
	particle = bins[i].particles;
	while (particle != NULL) {
	    move(*particle);
	    temp = particle;
	    particle = (*particle).next;
	    if (!isCorrectlyBinned((*temp).x,(*temp).y,bins[i].i,bins[i].j,bin_size)) {
		(*temp).next = rebinQueue;
		rebinQueue = temp;
	    } else {
		(*temp).next = stayQueue;
		stayQueue = temp;
	    }
	}
	bins[i].particles = stayQueue;
	stayQueue = NULL;
    }

    // rebin particles
    while (rebinQueue != NULL) {
	j = (int) ((*rebinQueue).x/bin_size);
	i = (int) ((*rebinQueue).y/bin_size);
	temp = rebinQueue;
	rebinQueue = (*rebinQueue).next;
	(*temp).next = bins[i*N+j].particles;
	bins[i*N+j].particles = temp;
    }
}

// debugging function
void check_bins(bin_t *bins, int N)
{
    particle_t *particle;
    particle_t *temp;
    int N2 = N*N;
    int counter = 1;
    for (int i = 0; i < N2; i++) {
	printf("%d,%d\n",bins[i].i,bins[i].j);
	particle = bins[i].particles;
	while (particle != NULL) {
	    printf("x=%f,y=%f\n",(*particle).x,(*particle).y);
	    temp = (*particle).next;
	    particle = temp;
	    counter++;
	    
        }
    }
    printf("there are %d particles in the bins\n",counter);
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
