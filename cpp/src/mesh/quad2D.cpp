#include "quad2D.h"

Quad2D::Quad2D(int id_,
               Elas2D* material_,
               const std::array<Node*,4>& nodes_) {
  id = id_;
  material = material_;
  for (int i = 0; i < 4; i++)
    nodes[i] = nodes_[i];
}