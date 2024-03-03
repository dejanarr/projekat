#ifndef IP_H
#define IP_H
#define SC_INCLUDE_FX

#include <iostream>
#include <systemc>
#include <string>
#include <fstream>
#include <deque>
#include <vector>
#include <array>
#include <algorithm>
#include "types.hpp"
#include "tlm_utils/tlm_quantumkeeper.h"

using namespace std;
using namespace sc_core;

SC_MODULE(Ip)
{

    public:
        SC_HAS_PROCESS(Ip);
        Ip(sc_module_name name);
        tlm_utils::simple_target_socket<Ip> s_ip_t0;
        
        sc_out<sc_logic> p_out;
        
    protected:
        void b_transport0(pl_t&, sc_time&);
        void proc();
        
        void AddSample(num_i r, num_i c, num_f rpos, num_f cpos, num_f rx, num_f cx, num_i step);
        void PlaceInIndex(num_f mag1, num_i ori1, num_f mag2, num_i ori2, num_f rx, num_f cx);
        
        vector <num_i> r;
        vector <num_i> c;
        vector <num_f> rpos;
        vector <num_f> cpos;
        vector <num_f> rx;
        vector <num_f> cx;
        vector <num_i> step;
        //_Pixels
        //_cose
        //_sine
        //_index
        
};

        
#endif
