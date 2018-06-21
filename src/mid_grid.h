
#ifndef __mid_grid__
#define __mid_grid__

#include "ising_solver.h"
#include "mid.h"
#include <vector>

class MidWithGrid : public Mid {
private:
  const int grid_size;
  // isingの各頂点が, 何ステップ目にどの頂点を訪れるものに対応するか
  std::vector<std::pair<int, int>> ising_node;
public:
  MidWithGrid(const Problem& prob, int grid_size);
  ~MidWithGrid();
  CostFunction getCostFunction();
  Answer getAnswerFromSpin(const std::vector<int>& spin) const;
};

#endif