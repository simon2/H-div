#include <cilk/cilk.h>
#include <cilk/reducer_list.h>
#include <iostream>

int main(){
  cilk::reducer<cilk::op_list_append<char> > result;
  cilk_for(std::size_t i = 'A'; i<'z'+1; ++i){
    result -> push_back((char)i);
  }

  std::cout << "String = ";
  std::list<char> r;
  r = result.get_value();
  for(std::list<char>::iterator i = r.begin();i != r.end(); ++i){
    std::cout << *i;
  }
  std::cout << std::endl;
}
