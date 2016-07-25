/*                             syscalls_ansi.c
                              -----------------
*/

/* Date:   19-Mar-2010
   Author: Bourdin.KIS (Bourdin@KIS.Uni-Freiburg.de)
   Description:
 ANSI C and standard library callable function wrappers for use in Fortran.
 Written to compensate for inadequatenesses in the Fortran95/2003 standards.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "headers_c.h"

/* ---------------------------------------------------------------------- */

void FTNIZE(file_size_c)
     (char *filename, FINT *bytes)
/* Determines the size of a file.
   Returns:
   * positive integer containing the file size of a given file
   * -2 if the file could not be found or opened
   * -1 if retrieving the file size failed
*/
{
  struct stat fileStat;
  int file = -1;

  *bytes = -2;
  file = open (filename, O_RDONLY);
  if(file == -1) return;

  *bytes = -1;
  if(fstat (file, &fileStat) < 0) { close (file); return; }
  close (file);

  *bytes = fileStat.st_size;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(write_binary_file_c)
     (char *filename, FINT *bytes, char *buffer, FINT *result)
/* Writes a given buffer to a binary file.
   Returns:
   * positive integer containing the number of written bytes
   * -2 if the file could not be opened
   * -1 if writing the buffer failed
*/
{
  int file = -1;
  int written = 0;

  *result = -2;
  file = open (filename, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);

  if(file == -1) return;

  *result = -1;

  written = (int) write (file, buffer, (size_t) *bytes);
  close (file);
  if (written != *bytes) return;
  *result = written;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_pid_c)
     (FINT *pid)
/* Determines the PID of the current process.
   Returns:
   * integer containing the PID of the current process
   * -1 if retrieving the PID failed
*/
{
  pid_t result;

  *pid = -1;
  result = getpid ();
  if (result) *pid = (int) result;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_env_var_c)
     (char *name, char *value)
/* Gets the content of an environment variable.
   Returns:
   * string containing the content of the environment variable, if available
   * empty string, if retrieving the environment variable failed
*/
{
  char *env_var;

  env_var = getenv (name);
  if (env_var) strncpy (value, env_var, strlen (env_var));
}

/* ---------------------------------------------------------------------- */

void FTNIZE(directory_exists_c)
     (char *path, FINT *exists)
/* Checks for existence of a directory.
   Returns:
   * 1, if 'path' points to a directory
   * -1, on error
   * 0, otherwise
*/
{
  struct stat result;

  *exists = stat (path, &result);
  if (S_ISDIR (result.st_mode)) *exists = 1;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(is_nan_c)
     (REAL *value, FINT *result)
/* Determine if value is not a number.
   Returns:
   * 1, if value is not a number
   * 0, if value is a number
   * -1 on failure (value is neither float or double)
*/
{
  *result = -1;

  if (sizeof (*value) == sizeof (double)) *result = isnan ((double) *value);
  /*
    isnanf() is sometimes not available
    if (sizeof (*value) == sizeof (float)) *result = isnanf ((float) *value);
  */
  if (sizeof (*value) == sizeof (float)) *result = !(*value == *value);
}

/* ---------------------------------------------------------------------- */

void FTNIZE(system_c) (char *command)
/* Date:   04-Nov-2011
   Author: MR (matthias.rheinhardt@helsinki.fi)
   Description: function wrapper for ANSI C function system.
*/
{
  int res=system(command);
  if (res == -1) return; // some error handling is missing here [Bourdin.KIS]
}

/* ---------------------------------------------------------------------- */

