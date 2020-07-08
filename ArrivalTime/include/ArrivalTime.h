#ifdef __cplusplus
extern "C" {
#endif
#ifndef ARRIVAL_TIME_H
#define ARRIVAL_TIME_H

typedef unsigned int BOOL;
#define TRUE 1
#define FALSE 0

#define STA_WIN_SIZE 4 
#define LTA_WIN_SIZE 32
#define STA_WIN_BIT 0
#define LTA_WIN_BIT 2
//#define RATIO_THRESHOLD 2.0
// #define BUFLEN 10000
#define CHS 3
#define DOWN_SAMPLING_FAC 10

#ifdef LINUX
typedef float dtype;
#else
typedef unsigned long dtype;
#endif

#define FRONTOFFSET 3
#define BACKOFFSET 7
#define PSSEPARATEWIN 1
#define MAXSAMPLINGRATE 500
#define MAXFILESECONDS 60*6*60
#define MAXNUMEVENTS 200 

typedef struct {
    int stationID;
    int samplingRate; // station or event ID, sampling rate per second
    int eventSeqNum;
    double x;
    double y;
    double z;
    double pTime;
    double sTime;
    double amplitude;
} EventInfo;

typedef struct {
    int samplingRate;
    unsigned int size;
    dtype *buffer0;
    dtype *buffer1;
    dtype *buffer2;
} DataWindow;

typedef struct {
    unsigned int size;
    unsigned int channelNum;
    int samplingRate;
    int stationID;
    double timingOffset;
    double x;
    double y;
    double z;
    //dtype *buffer[CHS];  
    dtype *buffer;  
} BufferInput;

typedef struct {
    int year;
    int month;
    int dayOfMonth;
} Date;

typedef struct {
    double minX;
    double minY;
    double minZ;
    double minLon;
    double minLat;
    double minSphrZ;
    double maxX;
    double maxY;
    double maxZ;
    double maxLon;
    double maxLat;
    double maxSphrZ;
} CoorHeader;

typedef struct {
    int eventIdx;
    double eventX;
    double eventY;
    double eventZ;
    double eventT;
} LocationInfo;

float mean(float *samples, int length);
float variance(float *samples, int length);
float variance_zmean(float *samples, int length);
BOOL detectEvent(DataWindow *dataWindow, double RATIO_THRESHOLD);
void separatePS(DataWindow *dataWindow, int separate_win_len, int channel, dtype *pbuffer, dtype *sbuffer);
int picking(dtype *buffer, int size, int pre_len);
void getArrivalTime(DataWindow *dataWindow, int separate_win_len, int channel, dtype *pbuffer, dtype *sbuffer, EventInfo *arrivalTime);
dtype maxap(dtype *buffer, int size);
void dispEventInfo(EventInfo event_info, int event_counter);
int dataProcess(BufferInput *bufInput, EventInfo *event_info, int id, double *streamStartTime, dtype *streamBuffer, double threshold);
//int dataProcess(int extraoffset, BufferInput *bufInput, EventInfo *event_info);

//void getArrivalTime(DataWindow *dataWindow, EventInfo *arrivalTime);

Date DayOfYearToDate(int year, int dayOfYear);

double ToEpochTime(int year, int month, int day, int hour, int minute, int second, int milliSecond, int dstOn);
void itoa(char *buf, int base, int d);
void WriteEventTxt(char *arrivalTimeFile, EventInfo *event, int eventCounter);

void WriteEventBin(char *arrivalTimeFile, EventInfo *event, int eventCounter, int id);

int DegToSec(double decimalDegree);

double SecToDeg(int totalSeconds);

#endif

#ifdef __cplusplus
}
#endif
