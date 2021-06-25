// The @GLIST class implements the order-based core maintenance algorithm
// in the paper. There is another class called @LLIST which has slightly
// better performance but the class is not included here.
#ifndef CORE_GLIST_GLIST_H_
#define CORE_GLIST_GLIST_H_

#include "core.h"

#include "heap.h"
#include "treap.h"

namespace core {
class GLIST final: public CoreMaintenance {
 public:
  explicit GLIST(const int n);
  ~GLIST();

  void ComputeCore(const std::vector<std::vector<int>>& graph,
                   const bool init_idx,
                   std::vector<int>& core);
  void Insert(const int v1, const int v2,
              std::vector<std::vector<int>>& graph,
              std::vector<int>& core);
  void Remove(const int v1, const int v2,
              std::vector<std::vector<int>>& graph,
              std::vector<int>& core);
  int CalcFollow(const int v, const int ___n,
              std::vector<std::vector<int>> graph,
              std::vector<int> core);     
  void Check(const std::vector<std::vector<int>>& graph,
             const std::vector<int>& core) const;

 private:
  struct ListNode {
    int rem;
    int ext;
    int prev;
    int next;
  };

  void Keep(const std::vector<std::vector<int>>& graph,
            const int v, const int K,
            const std::vector<int>& core,
            int& list_t, std::vector<int>& swap);
  void PropagateDismissal(const std::vector<std::vector<int>>& graph,
                          const int K, const int v,
                          std::vector<int>& core,
                          std::vector<int>& to_be_clear,
                          std::vector<int>& changed);
  void __PropagateDismissal(const std::vector<std::vector<int>>& graph,
                          const int K, const int v,
                          std::vector<int>& core,
                          std::vector<int>& to_be_clear,
                          std::vector<int>& changed);
  int GetRank(const int v) {
    if (0 == rank_[v]) {
      rank_[v] = tree_.Rank(v);
      garbage_.push_back(v);
    }
    return rank_[v];
  }
  int __GetRank(const int v) {
    if (0 == __rank_[v]) {
      __rank_[v] = __tree_.Rank(v);
      __garbage_.push_back(v);
    }
    return __rank_[v];
  }

  int n_,__n_;
  std::vector<int> head_,__head_;
  std::vector<int> tail_,__tail_;
  std::vector<ListNode> node_,__node_;
  std::vector<int> mcd_,__mcd_;
  std::vector<int> deg_,__deg_;
  std::vector<int> rank_,__rank_;
  std::vector<int> root_,__root_;
  std::vector<bool> evicted_,__evicted_;
  std::vector<bool> visited_,__visited_;
  gadget::Treap tree_,__tree_;
  gadget::MinHeap heap_,__heap_;
  std::vector<int> garbage_,__garbage_;
};
}  // namespace core

#endif
