#ifndef QUAD2D_H
#define QUAD2D_H

#include <iostream>
#include <vector>
#include <map>
#include "node.h"
#include "elas2D.h"

class Quad2D {
public:
  Quad2D(int id_,
         Elas2D* material_,
         const std::array<Node*,4>& nodes_);
  ~Quad2D();
public:
  int id;
  Elas2D* material;
  std::array<Node*,4> nodes;
  int intRule;
};

#endif // QUAD2D_H