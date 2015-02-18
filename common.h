#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

#define LEFT 1
#define RIGHT 2
#define UP 4
#define DOWN 8

//
// particle data structure
//
struct particle_t
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  unsigned int grid;
  unsigned int boundary;
  struct particle_t *next;
};

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void organize_particles( int n, particle_t *p, particle_t **first, unsigned int num_blocks );
particle_t** get_neighboring_particles( particle_t *first, unsigned int flag, unsigned int &n );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


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
