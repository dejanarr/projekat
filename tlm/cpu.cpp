#ifndef CPU_C
#define CPU_C
#include "cpu.hpp"



Cpu(sc_module_name name, char* image_file_name, char* labels_file_name):sc_module(name), image_file_name(image_file_name, labels_file_name(labels_file_name)
{
    SC_THREAD(software);
    p_port0.bind(sig0);
    sig0 = sc_dt::SC_LOGIC_0;
    labels.clear();
	
    cout << "Cpu constucted" << endl;
}

void Cpu:software()
{
//Ovo ce morati drugacije, za sada nek stoji tu
    int samplingStep = 2; // Initial sampling step (default 2)
    int octaves = 4; // Number of analysed octaves (default 4)
    double thres = 4.0; // Blob response treshold
    bool doubleImageSize = false; // Set this flag "true" to double the image size
    int initLobe = 3; // Initial lobe size, default 3 and 5 (with double image size)
    int indexSize = 4; // Spatial size of the descriptor window (default 4)
    struct timezone tz; struct timeval tim1, tim2; // Variables for the timing measure
    unsigned char *buf; //za prebacivanje podataka iz memorije u IP
    
    ImLoad ImageLoader;
    int arg = 0;
    string fn = "../data/out.surf";
    Image *im=NULL;
    while (++arg < argc) {
      if (! strcmp(argv[arg], "-i"))
        im = ImageLoader.readImage(argv[++arg]);
      if (! strcmp(argv[arg], "-o"))
        fn = argv[++arg];
    }
    
      // Start measuring the time
    gettimeofday(&tim1, &tz);

    // Create the integral image
    Image iimage(im, doubleImageSize);

    // Start finding the SURF points
      cout << "Finding SURFs...\n";

    // These are the interest points
    vector< Ipoint > ipts;
    ipts.reserve(300);

    // Extract interest points with Fast-Hessian
    FastHessian fh(&iimage, /* pointer to integral image */
                   ipts,
                   thres, /* blob response threshold */
                   doubleImageSize, /* double image size flag */
                   initLobe * 3 /* 3 times lobe size equals the mask size */,
                   samplingStep, /* subsample the blob response map */
                   octaves /* number of octaves to be analysed */);


    fh.getInterestPoints();

    // Initialise the SURF descriptor
    Surf des(&iimage, /* pointer to integral image */
             doubleImageSize, /* double image size flag */
             indexSize /* square size of the descriptor window (default 4x4)*/);

    // Get the length of the descriptor vector resulting from the parameters
    VLength = des.getVectLength();

    // Compute the orientation and the descriptor for every interest point
    for (unsigned n=0; n<ipts.size(); n++){
      // set the current interest point
      des.setIpoint(&ipts[n]);
      // assign reproducible orientation
      des.assignOrientation();
      // make the SURF descriptor
      des.makeDescriptor(); //kad pozovem ovu funkciju treba da procitam r odatle i da saljem u hardver
    }

    // stop measuring the time, we're all done
    gettimeofday(&tim2, &tz);

    // save the interest points in the output file
    saveIpoints(fn, ipts);

    // print some nice information on the command prompt
      cout << "Detection time: " <<
        (double)tim2.tv_sec + ((double)tim2.tv_usec)*1e-6 -
        (double)tim1.tv_sec - ((double)tim1.tv_usec)*1e-6 << endl;

/*********************************************************************************/
    
    pl_t pl;
    sc_time offset=SC_ZERO_TIME;

    #ifdef QUANTUM
    tlm_utils::tlm_quantumkeeper qk;
    qk.reset();
    #endif
	
    #ifdef QUANTUM
    qk.inc(sc_time(CLK_PERIOD, SC_NS));
    offset = qk.get_local_time();
    qk.set_and_sync(offset);
    #else
    offset += sc_time(CLK_PERIOD, SC_NS);
    #endif 

    ddr.clear();
    
    //***********************************************************************************
    //Proces slanja r iz softvera u hardver
    
    
    //slanje komande da se r smesti
    
    command = 0b00000001;
    buf = (unsigned char*)&command;
    
    pl.set_address(0x80000000);
    pl.set_data_length(sizeof(r));
    pl.set_command(TLM_WRITE_COMMAND);
    pl.set_data_ptr(buf);
    s_cp_i0->b_transport(pl, offset);
    assert(pl.get_response_status() == TLM_OK_RESPONSE);
    qk.set_and_sync(offset);		

    #ifdef QUANTUM
    qk.inc(sc_time(CLK_PERIOD, SC_NS));
    offset = qk.get_local_time();
    qk.set_and_sync(offset);
    #else
    offset += sc_time(CLK_PERIOD, SC_NS);
    #endif;
    
    /*ddr = ;
    buf=(unsigned char*)&ddr[0];

    pl.set_address(0);
    pl.set_data_length(SIZE_LETTER_DATA);
    pl.set_command(TLM_WRITE_COMMAND);
    pl.set_data_ptr(buf);
    s_cp_i1->b_transport(pl, offset); // po ovome znam gde se salje
    assert(pl.get_response_status() == TLM_OK_RESPONSE);
    qk.set_and_sync(offset);		

    ddr.clear();*/

}

//Morace da se izmeni, za sada nek stoji tu
void Cpu::saveIpoints(string sFileName, const vector< Ipoint >& ipts)
{
    ofstream ipfile(sFileName.c_str());
    if( !ipfile ) {
      cerr << "ERROR in loadIpoints(): "
           << "Couldn't open file '" << sFileName << "'!" << endl;
      return;
    }
  
  
    double sc;
    unsigned count = ipts.size();

    // Write the file header
    ipfile << VLength + 1 << endl << count << endl;

    for (unsigned n=0; n<ipts.size(); n++){
      // circular regions with diameter 5 x scale
      sc = 2.5 * ipts[n].scale; sc*=sc;
      ipfile  << ipts[n].x /* x-location of the interest point */
              << " " << ipts[n].y /* y-location of the interest point */
              << " " << 1.0/sc /* 1/r^2 */
              << " " << 0.0     //(*ipts)[n]->strength /* 0.0 */
              << " " << 1.0/sc; /* 1/r^2 */

      // Here should come the sign of the Laplacian. This is still an open issue
      // that will be fixed in the next version. For the matching, just ignore it
      // at the moment.
      ipfile << " " << 0.0; //(*ipts)[n]->laplace;

      // Here comes the descriptor
      for (int i = 0; i < VLength; i++) {
        ipfile << " " << ipts[n].ivec[i];
      }
      ipfile << endl;
    }

    // Write message to terminal.
      cout << count << " interest points found" << endl;
}

