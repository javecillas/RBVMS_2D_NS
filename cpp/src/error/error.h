// ---- error.h ---- Fri Nov 5 09:06:22 CST 2004
// (c) Copyright: C. Armando Duarte, 2002-2004. All rights reserved
// See README.CopyRight file for further details.
//

#ifndef ERRORH
#define ERRORH

#include <iostream>
#include <string>

/**
 *  Prints a message with the file and line where the error occured
 *
 *  @param line input: line of the error
 *  @param file input: file of the error
 *  @param out input: type of output for message. Default is std::cerr
 */
void errLoc( const int line, const std::string &file, 
             std::ostream &out =std::cerr);

/**
 *  DEPRECATED: Old error stop function that shows file and line where
 *              error occured and then exit(EXIT_FAILURE) program
 *
 *  @param message input: message to be shown before exiting program
 *  @param out input: type of output for message. Default is std::cerr
 */
void errStop( const std::string &message, std::ostream &out = std::cerr);

/**
 *  In debugger, throwing an exception makes the debugger stop exactly where the
 *  ERRCHK or ERREXIT happened showing all the stack.
 *  This functions prints the file and line of the error and then
 *  throws an exception.
 *
 *  @param message input: message to be shown before throwing exception
 *  @param out input: type of output for message. Default is std::cerr
 */
void DebugStop( const std::string &message, std::ostream &out = std::cerr);

// *******************************************************************
// *******************************************************************


/**
 Check if a file exists.
 NS: Even though this function is not an "error" function per say
 I'm putting it here so everyone can see it

 @param name name of the file to check if exists
 @return if file exists
 */
bool file_exists (const std::string& name);

// New error function
// Additionally to the file and error line, it throws an exception.
// This is particularly good for debuggers because it shows all
// the stack of functions called before the bug.
// In gcc and LLVM, the error message is the same as exit(EXIT_FAILURE)

#define ERREXIT(str)  do{errLoc(__LINE__,__FILE__);DebugStop(str);}while(0)

#define ERRCHK(ex,str)   do { if(ex) { errLoc(__LINE__,__FILE__); \
                                       DebugStop(str);   }  } while(0)


// Old error function.
// it only prints the file and line of the error

/*
#define ERREXIT(str)  do{errLoc(__LINE__,__FILE__);errStop(str);}while(0)

#define ERRCHK(ex,str)   do { if(ex) { errLoc(__LINE__,__FILE__); \
                                       errStop(str);   }  } while(0)
*/




// *******************************************************************
// *******************************************************************

#endif
