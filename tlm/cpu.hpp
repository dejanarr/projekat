#ifndef CPU_H
#define CPU_H
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
#include "addresses.hpp"
#include "tlm_utils/tlm_quantumkeeper.h"
#include "../spec/imload.h"
#include "../spec/image.h"
#include "../spec/fasthessian.h"
#include "../spec/surf.h"

using namespace std;
using namespace sc_core;

SC_MODULE(Cpu)
{

	public:
		SC_HAS_PROCESS(Cpu);
		Cpu(sc_module_name name, char* image_file_name);
		tlm_utils::simple_initiator_socket<Cpu> s_cp_i0;
		tlm_utils::simple_initiator_socket<Cpu> s_cp_i1;
		sc_out<sc_dt::sc_logic> p_port0;
		sc_signal<sc_dt::sc_logic> sig0;	

	protected:
		void software();
};


#endif
