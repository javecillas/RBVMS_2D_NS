#include <fstream>
#include "readsolmech.h"

using namespace boost;

// **************************************************************
// **************************************************************

ReadSolMech::ReadSolMech(const std::string& fileIn) {
  // Test that input file opens properly
  std::ifstream infile(fileIn);
  if (!infile){
    std::cerr << "Unable open input file \"" 
              << fileIn << "\"" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // Start reading the input file
  std::string lineString, section;
  const std::string mysep = ",= ";
  while(!infile.eof()){
    if (section == "*NODE" || section == "*Node") {
    } else {
      getline(infile, lineString);
      section = getFirstWordWithSeparator(mysep,lineString);
    }
  }
}

// **************************************************************
// **************************************************************

ReadSolMech::~ReadSolMech() {
}

// **************************************************************
// **************************************************************

PhTokenizer 
ReadSolMech::getTokenizer(const std::string& mysep,
                          const std::string& lineString) const {
  boost::char_separator<char> sep(mysep.c_str());
  PhTokenizer inputTokens(lineString, sep);
  PhTokenizer::iterator tok_iterz;
  inputTokens.assign(lineString);
  return inputTokens;
}

// **************************************************************
// **************************************************************

const std::string 
ReadSolMech::getFirstWordWithSeparator(const std::string &mysep, 
                                       const std::string &lineString) const {
  PhTokenizer inputTokens = getTokenizer(mysep,lineString);
  PhTokenizer::iterator tok_iterz = inputTokens.begin();
  std::string returnstring;
  if (tok_iterz == inputTokens.end()) {
    returnstring = "nothing";
  } else {
    returnstring = *tok_iterz;
  }
  return returnstring;
}