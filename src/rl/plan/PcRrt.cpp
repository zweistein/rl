#include <iostream>

#include "PcRrt.h"

namespace rl
{
  namespace plan
  {
    ::std::string PcRrt::getName() const 
    {
      return "PCRRT";
    }

    bool PcRrt::solve() 
    {
      ::std::cout << "PcRrt solve!" << ::std::endl;
      return Rrt::solve();
    }
  }
}