#include <sys/time.h>

double wclock_()
{
  struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
