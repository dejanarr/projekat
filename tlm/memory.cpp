#ifndef MEM_C
#define MEM_C
#include "memory.hpp"

Mem::Mem(sc_module_name name):sc_module(name)
{
    s_mem_t0.register_b_transport(this, &Mem::b_transport);
    s_mem_t1.register_b_transport(this, &Mem::b_transport);
	
    cout<<"Memory constructed"<<endl;
    ddr.reserve(1333*2000*3);
} 

void Mem::b_transport(pl_t& pl, sc_time& offset)
{
    tlm_command cmd    = pl.get_command();
    sc_dt::uint64 adr  = pl.get_address();
    unsigned char *buf = pl.get_data_ptr();
    unsigned int len   = pl.get_data_length();
        
    switch(cmd)
    {
        case TLM_WRITE_COMMAND:
            for(unsigned int i=0; i<len; i++)
            {       
                ddr[adr+i]=((num_f*)buf)[i];
            }
            pl.set_response_status(TLM_OK_RESPONSE);
            break;

        case TLM_READ_COMMAND:
            buf = (unsigned char*)&ddr[adr];
            pl.set_data_ptr(buf);
            pl.set_response_status(TLM_OK_RESPONSE);
            break;
        
        default:
            pl.set_response_status( TLM_COMMAND_ERROR_RESPONSE );
    }
    
    offset += sc_time(10, SC_NS);
}

#endif // MEM_C





/*#include <cstring> // Za memcpy

void Mem::b_transport(pl_t& pl, sc_time& offset) {
    tlm_command cmd = pl.get_command();
    sc_dt::uint64 adr = pl.get_address();
    unsigned char* buf = pl.get_data_ptr();
    unsigned int len = pl.get_data_length();

    switch (cmd) {
        case TLM_WRITE_COMMAND: {
            // Pretpostavimo da je len broj bajtova, a ne broj elemenata.
            unsigned int num_elements = len / sizeof(num_f); // Ili sizeof(num_i) zavisno od konteksta
            for (unsigned int i = 0; i < num_elements; i++) {
                // Pretpostavimo da ddr ima dovoljno prostora za skladištenje
                // Ovde treba biti oprezan zbog prelivanja adresa
                memcpy(&ddr[adr / sizeof(num_f) + i], buf + i * sizeof(num_f), sizeof(num_f));
            }
            pl.set_response_status(TLM_OK_RESPONSE);
            break;
        }
        case TLM_READ_COMMAND: {
            // Slično kao za pisanje, ali obrnuto
            unsigned int num_elements = len / sizeof(num_f); // Ili sizeof(num_i) zavisno od konteksta
            for (unsigned int i = 0; i < num_elements; i++) {
                memcpy(buf + i * sizeof(num_f), &ddr[adr / sizeof(num_f) + i], sizeof(num_f));
            }
            pl.set_data_ptr(buf);
            pl.set_response_status(TLM_OK_RESPONSE);
            break;
        }
        default:
            pl.set_response_status(TLM_COMMAND_ERROR_RESPONSE);
    }
    
    offset += sc_time(10, SC_NS);
}

Napomene:
Ovaj pristup pretpostavlja da su svi podaci koji se čitaju ili pišu istog tipa (num_f ili num_i). Ako planirate da imate memorijsku mapu koja sadrži oba tipa, možda će biti potrebno da implementirate dodatnu logiku za upravljanje tipovima podataka.
Takođe, bitno je napomenuti da memcpy radi na nivou bajtova, što znači da treba pažljivo obratiti pažnju na veličinu i poravnanje podataka kako bi se izbegli problemi sa pristupom memoriji.
Korišćenje / sizeof(num_f) (ili num_i) za izračunavanje indeksa u ddr nizu pretpostavlja da je svaki element u ddr nizu tačno veličine num_f (ili num_i). Ovo može zahtevati dodatnu logiku ako ddr sadrži mešovite tipove podataka.


*/

//BOLJE

/*#include <variant>
#include <vector>

// Definišite varijantni tip koji može da sadrži int i double
using IntOrDouble = std::variant<int, double>;

class MixedMemory {
public:
    std::vector<IntOrDouble> memory;

    void write(size_t address, IntOrDouble value) {
        // Proverite granice i upišite vrednost
        if (address < memory.size()) {
            memory[address] = value;
        }
    }

    IntOrDouble read(size_t address) {
        // Proverite granice i pročitajte vrednost
        if (address < memory.size()) {
            return memory[address];
        }
        return {}; // Vraća prazan std::variant kao grešku
    }
};

