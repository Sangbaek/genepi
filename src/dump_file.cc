#include "cpp.h"

#include "hepevt_genepi.h"
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
    fprintf(ptr,"%5d %5d %5d %5d %5d %5d %5.1f %5d %5d %15.10f\n",
            3,trk.TarA, trk.TarZ, trk.Theli, trk.Bheli, 11, trk.Eb, trk.Type[1], 1, xsec);
    TRandom1 rndm_vertex;
    double vx=rndm_vertex.Uniform(0,1)*0.75;
    double vy=rndm_vertex.Uniform(0,1)*0.75;
    double vz=0;
    for(int i=0; i<trk.Ntracks; i++)
    {
      int index;
      if (i==2) index=1;
      else if (i==5) index=2;
      else if (i==6) index=3;
      else continue;
      fprintf(ptr,"%5d %5d %5d %5d %5d %5d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
          index, trk.Charge[i],1, trk.Type[i],0,0, trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i],sqrt(sqr(trk.E[i])-sqr(trk.P[i])),vx,vy,vz);
      }
  }

}
