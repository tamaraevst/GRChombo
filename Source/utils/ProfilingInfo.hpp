// MK: This class implements the possibility of basic profiling info for
// GRChombo.  It worked as a quick profiling hook for PAPI but not much thought
// has been put into this.

#ifndef PROFILINGINFO_HPP_
#define PROFILINGINFO_HPP_

class ProfilingInfo
{
  public:
    ProfilingInfo(){};
    virtual ~ProfilingInfo(){};

    virtual bool isActive() = 0;

    virtual int startCounters() = 0;
    virtual int readCounters() = 0;
    virtual int readShutdown() = 0;
    virtual int resetCounters() = 0;
    virtual int readTotalTimes() = 0;
};
#endif /* PROFILINGINFO_HPP_ */
