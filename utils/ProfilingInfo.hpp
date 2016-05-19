//MK: This class implements the possibility of basic profiling info for GRChombo.
//It works as a quick profiling hook but not much though has been put into this.

#ifndef _PROFILINGINFO_HPP
#define _PROFILINGINFO_HPP

class ProfilingInfo
{
   public:
      ProfilingInfo(){};
      ~ProfilingInfo();

      virtual bool isActive() = 0;

      virtual int startCounters() = 0;
      virtual int readCounters() = 0;
      virtual int readShutdown() = 0;
      virtual int resetCounters() = 0;
      virtual int readTotalTimes() = 0;
};
#endif //_PROFILINGINFO_HPP
