// ---- error.C ---- Fri Nov 5 09:06:15 CST 2004
// (c) Copyright: C. Armando Duarte, 2002-2004. All rights reserved
// See README.CopyRight file for further details.
//


#include <cstdio>
#include <cstdlib>
#include "error.h"
#include <unistd.h>

using std::endl;

void errLoc( const int line, const std::string &file, std::ostream &out ) {

  out << "In line " << line << " of file " << file << ": ";
}
// *******************************************************************
// *******************************************************************

void errStop( const std::string &message, std::ostream &out ) {

  out << message << endl;
  exit(EXIT_FAILURE);
}
// *******************************************************************
// *******************************************************************

void DebugStop( const std::string &message, std::ostream &out ) {
  
  out << message << endl;
  std::bad_exception bad;
  throw bad;
}


// *******************************************************************
// *******************************************************************

bool boost_error(char const * expr, char const * func, char const * file, 
                 long line)
{
  std::printf("%s(%ld): Assertion '%s' failed in function '%s'\n", file, 
              line, expr, func);
  return true; // fail w/ standard assert()
}
// *******************************************************************
// *******************************************************************

bool file_exists (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
