#ifndef MEM_H
#define MEM_H

#include <iostream>
#include <systemc>
#include <string>
#include <fstream>
#include <vector>
#include "types.hpp"
#include "adresses.hpp"
#include "tlm_utils/tlm_quantumkeeper.h"

using namespace std;
using namespace sc_core;

SC_MODULE(Mem)
{
public:
    //SC_HAS_PROCESS(Mem);
    Mem(sc_module_name name);

    tlm_utils::simple_target_socket<Mem> s_mem_t0; 
    tlm_utils::simple_target_socket<Mem> s_mem_t1;
    void b_transport(pl_t&, sc_time&);

protected:
    vector<num_f> ddr;		
};

#endif // MEM_H
