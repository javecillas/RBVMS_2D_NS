#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>

#include <array>

class Node {
public:
  // Constructor 
  Node(int id_, const std::array<double,2>& coor_);
  // Destructor
  ~Node();
public:
  // Node ID
  int id;
  // Node coordiante
  std::array<double,2> coor;
  // Node global DOF
  std::array<int,2> gDoF;
  // Flag
  bool isDoFSet;
  bool hasDirBC;
  double bcVal;
  int bcFlag;
};

#endif // NODE_H