// Ref: Data Clustering in C++ - Guojun Gan
#include<set>

class Cluster {
public:
  Cluster(int id): _id(id){
  }
  void add_item(int item) {
    _items.insert(item);
  }
  void remove_item(int item){
    _items.erase(item);
  }
  int size() const {
    return _items.size();
  }
  const std::set<int>& data() const {
    return _items;
  }
private:
  int _id;
  std::set<int> _items;
};
