#include "cpp.h"

#include "hepevt.h"
#include "track_vars.h"
#include "dump_file.h"

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
}
