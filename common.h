#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//

typedef struct particle_t
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  particle_t *next;
} particle_t;


typedef struct bin_t
{
  int i; // row index of block
  int j; // col index of block
  bin_t *north;
  bin_t *south;
  bin_t *east;
  bin_t *west;
  bin_t *ne;
  bin_t *se;
  bin_t *nw;
  bin_t *sw;
  particle_t *particles;
} bin_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
double get_size();
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );

void initialize_bins(bin_t *bins, int N);
void bin_particles(bin_t *bins, particle_t *particles, int N, int n, double bin_size);
void apply_force_bin(particle_t *self, particle_t *neighbor, bool clearForces, double*dmin, double *davg, int *navg);
void interact_binned_particles(bin_t *bins, int N, double *dmin, double *davg, int *navg);
bool isCorrectlyBinned(double x, double y, int i, int j, double bin_size);
void move_and_rebin(bin_t *bins, int N, double bin_size);
void check_bins(bin_t *bins, int N);
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
