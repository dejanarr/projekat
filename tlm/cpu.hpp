#ifndef CPU_H
#define CPU_H
#define SC_INCLUDE_FX

#include <iostream>
#include <sys/time.h>
#include <string.h>
#include <cmath>
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

using namespace std;
using namespace sc_core;

typedef sc_dt::sc_int<11> num_i;
typedef sc_dt::sc_fixed<48, 30, sc_dt::SC_TRN, sc_dt::SC_SAT> num_f;

SC_MODULE(Cpu)
{

	public:
		SC_HAS_PROCESS(Cpu);
		Cpu(sc_module_name name, char* image_file_name);
		tlm_utils::simple_initiator_socket<Cpu> s_cp_i0;
		tlm_utils::simple_initiator_socket<Cpu> s_cp_i1;
		sc_out<sc_dt::sc_logic> p_port0;
		sc_signal<sc_dt::sc_logic> sig0;
		
		
		Image *_iimage = nullptr;
		Ipoint *_current = nullptr;
		std::vector<std::vector<std::vector<num_f>>> _index;
		bool _doubleImage = false;
		num_i _VecLength = 0;
		num_i _IndexSize = 4;
		num_i _MagFactor = 0;
		num_i _OriSize = 0;
		num_i _width = 0, _height = 0;

		num_f _sine = 0.0, _cose = 1.0;
		std::vector<std::vector<num_f>> _Pixels;

		num_f _lookup1[83], _lookup2[40];
		
		int VLength; // Length of the descriptor vector
		
		double scale;
		int x;
		int y;	

	protected:
		void software();
		void createVector(double scale, double row, double col);
		void normalise();
		void createLookups();
		void initializeGlobals(Image *im, bool dbl, int insi);
		int getVectLength();
		void setIpoint(Ipoint* ipt);
		void assignOrientation();
		void makeDescriptor();
		void saveIpoints(string sFileName, const vector< Ipoint >& ipts);
	
};


#endif
