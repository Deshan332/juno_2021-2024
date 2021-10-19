#include <string>
#include <time.h>
using namespace std;
TString timestamp()
{
    time_t rawtime;
    struct tm *info;
    char run_time[80];
    time( &rawtime );
    info = localtime( &rawtime );
    strftime(run_time,80,"%F_%I:%M", info);
    return run_time;
}