#ifndef READSOLMECH_H
#define READSOLMECH_H

#include <iostream>
#include <string>
#include <boost/tokenizer.hpp>

typedef boost::tokenizer<boost::char_separator<char>> PhTokenizer;

class ReadSolMech {
public:
  ReadSolMech(const std::string& fileIn);
  ~ReadSolMech();

  PhTokenizer getTokenizer(const std::string& mysep,
                           const std::string& lineString) const;
  const std::string 
  getFirstWordWithSeparator(const std::string &mysep, 
                            const std::string &lineString) const;

private:

};

#endif // SOLMECH_H