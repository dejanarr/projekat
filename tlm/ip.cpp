#ifndef IP_C
#define IP_C
#include "ip.hpp"


Ip::Ip(sc_module_name name) : sc_module(name)
{
    SC_THREAD(proc);
    s_ip_t0.register_b_transport(this, &Ip::b_transport0);
    
    command = 0b00000000;
    
    /*...*/
    
    cout << "IP created" << endl;
}

void Ip::b_transport0(pl_t& pl, sc_time& offset)
{
    tlm_command cmd = pl.get_command();
    sc_dt::uint64 adr = pl.get_adress();
    const unsigned char* buf = pl.get_data_ptr();
    unsigned int len = pl.get_data_length();
    
    sc_dt::uint64 temp = adr - 0x80000000;
    
    switch (cmd)
    {
        case TLM_WRITE_COMMAND:
            switch (temp)
            {
                case 0x00000001:
                    command = int(*buf);
                    pl.set_response_status(TLM_OK_RESPONSE);
                    break;
                case 0x00000010:
                    //...//
                    pl.set_response_status(TLM_OK_RESPONSE);
                    break;
                case 0x00000100:
                    //...//
                    pl.set_response_status(TLM_OK_RESPONSE);
                    break;
                default:
                    pl.set_response_status(TLM_COMMAND_ERROR_RESPONSE);
                }
                break;
                
        case TLM_READ_COMMAND:
            break;
            
        default:
            pl.set_response_status(TLM_COMMAND_ERROR_RESPONSE);
    }
    offset += sc_time(10, SC_NS);
}

void Ip::proc()
{

}


//Nek stoji za sada ovako
void Ip::AddSample(num_i r, num_i c, num_f rpos,
                     num_f cpos, num_f rx, num_f cx, num_i step) {
  num_f weight;
  num_f dx, dy;
  
  if (r < 1+step  ||  r >= _height - 1-step  ||
      c < 1+step  ||  c >= _width - 1-step)
     return;
 
  weight = _lookup2[num_i(rpos * rpos + cpos * cpos)];
  
  num_f dxx, dyy;

  dxx = weight*get_wavelet2(_Pixels, c, r, step);
  dyy = weight*get_wavelet1(_Pixels, c, r, step);
  dx = _cose*dxx + _sine*dyy;
  dy = _sine*dxx - _cose*dyy;

 PlaceInIndex(dx, (dx<0?0:1), dy, (dy<0?2:3), rx, cx);
 
}

void Ip::PlaceInIndex(num_f mag1, num_i ori1, num_f mag2, num_i ori2, num_f rx, num_f cx) {

  num_i ri = std::max(0, std::min(static_cast<int>(_IndexSize - 1), static_cast<int>(rx)));
  num_i ci = std::max(0, std::min(static_cast<int>(_IndexSize - 1), static_cast<int>(cx)));
  
  num_f rfrac = rx - ri;
  num_f cfrac = cx - ci;
  
  rfrac = std::max(0.0f, std::min(float(rfrac), 1.0f));
  cfrac = std::max(0.0f, std::min(float(cfrac), 1.0f));
  
  num_f rweight1 = mag1 * (1.0 - rfrac);
  num_f rweight2 = mag2 * (1.0 - rfrac);
  num_f cweight1 = rweight1 * (1.0 - cfrac);
  num_f cweight2 = rweight2 * (1.0 - cfrac);
  

  if (ri >= 0 && ri < _IndexSize && ci >= 0 && ci < _IndexSize) {
    _index[ri][ci][ori1] += cweight1;
    _index[ri][ci][ori2] += cweight2;
  }

  if (ci + 1 < _IndexSize) {
    _index[ri][ci + 1][ori1] += rweight1 * cfrac;
    _index[ri][ci + 1][ori2] += rweight2 * cfrac;
  }

  if (ri + 1 < _IndexSize) {
    _index[ri + 1][ci][ori1] += mag1 * rfrac * (1.0 - cfrac);
    _index[ri + 1][ci][ori2] += mag2 * rfrac * (1.0 - cfrac);
  }
}


#endif // IP_C
