#include "cpp.h"

#include "hepevt.h"
#include "track_vars.h"
#include "dump_file.h"
#include "root.h"
#include "genepi.h"
//void dump_file(int mode, double xsec, FILE *ptr)
void dump_file(int mode, double xsec)
{
  if(mode == 0)
  {
    fprintf(ptr,"%d\n",hepevt.NHEP);
    for(int i=0; i<hepevt.NHEP; i++)
    {
      fprintf(ptr,"%5d %5d %5d %5d %15.8f %15.8f %15.8f %15.8f\n",
                hepevt.ISTHEP[i], hepevt.IDHEP[i], hepevt.JDAHEP[i][0], 
		hepevt.JDAHEP[i][1], hepevt.PHEP[i][0], hepevt.PHEP[i][1], 
		hepevt.PHEP[i][2],hepevt.PHEP[i][4]);
    }
  }
  else if(mode == 1)
  {
    fprintf(ptr,"%6d %6d %6d %8.4f %8.4f %8.4f %8.4f %8.4f %15.10f\n",
            trk.Ntracks, trk.TarA, trk.TarZ, trk.xbj, trk.y, trk.W2, 
	    trk.Q2, trk.nu, xsec);
    for(int i=0; i<trk.Ntracks; i++)
    {
      fprintf(ptr,"%5d %5d %5d %15.8f %15.8f %15.8f %15.8f\n",
          i, trk.Charge[i], trk.Type[i], trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i]);
    }
  }
  // compatable with gemc
  else if(mode == 2)
  {
    fprintf(ptr,"%6d %6d %6d %8.1f %8.1f %6d %8.3f %6d %6d %15.10f\n",
            3, trk.TarA, trk.TarZ, trk.Theli, trk.Bheli, 11, 
      trk.Eb, trk.Type[1],1, xsec);
    TRandom1 rndm_vertex;
    rndm_vertex.SetSeed(12345);
    double vx=rndm_vertex.Rndm()*0.75;
    double vy=rndm_vertex.Rndm()*0.75;
    double vz=0;
    for(int i=0; i<trk.Ntracks; i++)
    {
    if (i!=2 && i!=5 && i!=6 ) continue;
    fprintf(ptr,"%5d %5d %5d %5d %5d %5d %15.8f %15.8f %15.8f %15.8f %5d %15.8f %15.8f %15.8f\n",
        i, trk.Charge[i],1, trk.Type[i],0,0, trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i],sqrt(sqr(trk.E)-sqr(trk.P)),vx,vy,vz);
    }
  }

}
