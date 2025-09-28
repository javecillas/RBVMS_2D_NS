#include "node.h"

Node::Node(int id_, const std::array<double,2>& coor_) {
  id = id_;
  coor = coor_;
  //
  isDoFSet  = false;
  hasDirBC  = false;
}

Node::~Node() {
}

