#ifndef ELAS2D_H
#define ELAS2D_H

#include <iostream>

class Elas2D {
public:
  // Constructor and destructor
  Elas2D(int id_, double E_, double nu_);
  ~Elas2D();
  // Member functions
  double getE() { return E; };
public:
  // Member variables
  int id;
  double E;
  double nu;
};

#endif // ELAS2D_H