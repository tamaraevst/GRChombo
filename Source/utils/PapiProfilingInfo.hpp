// MK: This class implements the possibility of basic profiling info for
// GRChombo using PAPI  It works as a quick profiling tool but not much though has
// been put into this.
#ifndef _PAPIPROFILINGINFO_HPP
#define _PAPIPROFILINGINFO_HPP
#ifdef USE_PAPI
#include "papi.h"
#endif

using std::string;

class PapiProfilingInfo : public ProfilingInfo
{
  public:
    PapiProfilingInfo(int *a_events, const int a_numEvents);

    ~PapiProfilingInfo();

    // Some getters
    bool isActive();
    long long *getCounters();
    long long getLastRealTime();
    long long getLastProcTime();

    void shutdown();

    int startCounters();
    int readCounters();
    int readShutdown();
    int resetCounters();
    int readTotalTimes();

    bool noError(int errorCode, bool shutdownOnError = true);

  private:
    bool usingPapi;
    int numEvents;
    char *eventNames;
    long long lastRealTime;
    long long lastProcTime;
    long long initialRealTime;
    long long initialProcTime;
    float sumReadRealSeconds; // This gives the sum of all times that were read
                              // out, in seconds
    float sumReadProcSeconds;
    long long *counters;
    int eventSet;
};

PapiProfilingInfo::PapiProfilingInfo(int *a_events, int a_numEvents)
    : numEvents(a_numEvents), lastRealTime(0), lastProcTime(0), usingPapi(true),
      counters(new long long[a_numEvents]),
      eventNames(new char[a_numEvents * PAPI_MAX_STR_LEN]), eventSet(PAPI_NULL),
      sumReadProcSeconds(0), sumReadRealSeconds(0)
{
    // Throw a warning if Papi is already initialized
    if (PAPI_is_initialized())
        MayDay::Warning("Papi has already been initialized. Use of this class "
                        "is not safe in this case.");

    // Initialize Papi
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT)
    {
        printf("PAPI_library_init returned: %d. Papi version is: %d. \n",
               retval, PAPI_VER_CURRENT);
        shutdown();
    }
    else
    {
        for (int i = 0; i < a_numEvents; ++i)
        {
            // Fill in the name of the event into event name array
            PAPI_event_code_to_name(a_events[i],
                                    eventNames + i * PAPI_MAX_STR_LEN);
            // Check whether the event exists on the architecture
            if (PAPI_query_event(a_events[i]) != PAPI_OK)
            {
                printf("Event %c ist not available.",
                       (eventNames + i * PAPI_MAX_STR_LEN));
                shutdown();
            }
        }
        pout() << "Using papi..." << endl;
        // Add the events to the eventSet
        int retval = PAPI_create_eventset(&eventSet);
        if (noError(retval))
        {
            int retval = PAPI_add_events(eventSet, a_events, a_numEvents);
            noError(retval);
        }
    }
}

PapiProfilingInfo::~PapiProfilingInfo()
{
    // shutdown papi
    shutdown();
    delete[] counters, eventNames;
}

/*Some getters*/
bool PapiProfilingInfo::isActive() { return usingPapi; }

long long *PapiProfilingInfo::getCounters() { return counters; }

long long PapiProfilingInfo::getLastRealTime() { return lastRealTime; }

long long PapiProfilingInfo::getLastProcTime() { return lastProcTime; }

// TODO: Add setters - one might want them

void PapiProfilingInfo::shutdown()
{
    usingPapi = false;
    PAPI_shutdown();
}

int PapiProfilingInfo::startCounters()
{
    // Return error if Papi isn't active anymore
    if (!isActive())
        return -1;

    int retval = PAPI_start(eventSet);
    if (!noError(retval))
        return -2;

    // Just before we start to run initialize the timers to the current time:
    initialRealTime = PAPI_get_real_usec();
    initialProcTime = PAPI_get_virt_usec();
    lastRealTime = initialRealTime;
    lastProcTime = initialProcTime;

    return 1;
}

int PapiProfilingInfo::readCounters()
{
    // Return error if Papi isn't active anymore
    if (!isActive())
        return -1;

    /* Read elapsed real and process times  */
    long long rt = PAPI_get_real_usec();
    long long pt = PAPI_get_virt_usec();
    float real_time = (static_cast<float>(rt - lastRealTime)) * .000001;
    float proc_time = (static_cast<float>(pt - lastProcTime)) * .000001;
    sumReadRealSeconds += real_time;
    sumReadProcSeconds += proc_time;

    // Output the timing
    pout() << "Real time   " << real_time << "   Proc_time " << proc_time;

    int retval = PAPI_read(eventSet, counters);
    // Collect the counter values
    if (!noError(retval))
        return -2;
    else
    {
        for (int i = 0; i < numEvents; ++i)
        {
            pout() << "   " << (eventNames + i * PAPI_MAX_STR_LEN) << "   "
                   << counters[i];
        }
        pout() << endl;
    }
    return 1;
}

int PapiProfilingInfo::resetCounters()
{
    if (!isActive())
        return -1; // Return error if Papi isn't active anymore

    int retval = PAPI_reset(eventSet);
    lastRealTime = PAPI_get_real_usec();
    lastProcTime = PAPI_get_virt_usec();
    return noError(retval, false);
}

int PapiProfilingInfo::readShutdown()
{
    if (!isActive())
        return -1; // Return error if Papi isn't active anymore
    pout() << " --- Final PAPI counter output --- " << endl;

    /* Read elapsed real and process times  */
    long long rt = PAPI_get_real_usec();
    long long pt = PAPI_get_virt_usec();
    float real_time = (static_cast<float>(rt - lastRealTime)) * .000001;
    float proc_time = (static_cast<float>(pt - lastProcTime)) * .000001;
    float totalRealTime = (static_cast<float>(rt - initialRealTime)) * .000001;
    float totaProcTime = (static_cast<float>(pt - initialProcTime)) * .000001;

    // Output the differential timing
    pout() << "Real time   " << real_time << "   Proc_time " << proc_time;

    int retval = PAPI_stop(eventSet, counters);
    // Collect the counter values
    if (!noError(retval))
        return -2;
    else
    {
        for (int i = 0; i < numEvents; ++i)
        {
            pout() << "   " << (eventNames + i * PAPI_MAX_STR_LEN) << "   "
                   << counters[i];
        }
        pout() << endl;
    }

    readTotalTimes();

    shutdown();
    return 1;
}

int PapiProfilingInfo::readTotalTimes()
{
    if (!isActive())
        return -1; // Return error if Papi isn't active anymore

    /* Read elapsed real and process times  */
    long long rt = PAPI_get_real_usec();
    long long pt = PAPI_get_virt_usec();
    float totalRealTime = (static_cast<float>(rt - initialRealTime)) * .000001;
    float totalProcTime = (static_cast<float>(pt - initialProcTime)) * .000001;
    float fracTime[2] = {sumReadRealSeconds / totalRealTime,
                         sumReadProcSeconds / totalProcTime};

    pout() << "Total real time   " << totalRealTime << "   Total proc time "
           << totalProcTime << endl;
    pout() << "Sum of read out values: real time " << sumReadRealSeconds
           << " proc time " << sumReadProcSeconds << endl;
    pout() << "As fraction of total time real time: " << fracTime[0]
           << " proc time " << fracTime[1] << endl;

#ifdef CH_MPI
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // If we are using MPI also average the fraction above over the ranks
    float procSums[2];
    MPI_Reduce(&fracTime, &procSums, 2, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
        pout() << "The above averaged over all ranks: real time "
               << procSums[0] / size << " proc time " << procSums[1] / size
               << endl;
#endif
    return 1;
}

bool PapiProfilingInfo::noError(int errorCode, bool shutdownOnError)
{
    if (errorCode < PAPI_OK)
    {
        printf("PAPI internal error %d: %s \n", errorCode,
               PAPI_strerror(errorCode));
        if (shutdownOnError)
        {
            printf("Shutting down PAPI... \n");
            shutdown();
        }
        return false;
    }
    return true;
}
#endif //_PAPIPROFILINGINFO_HPP
