#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <math.h>
#include <unistd.h>
#include <utility>
#include "ArrivalTime.h"
#include "eigen.h"
#include "matrix.h"

#ifdef LINUX
#include <math.h>
#else
#include "QmathLib.h"
#endif
 
using namespace std;
#ifdef LINUX
#define LOG(x) log(x)
#else
#define GLOBAL_Q 12
#define LOG(x) _QtoF(_Qlog(_Q(sigma_2)))
#endif

char outputTxtFile[] = "a.csv";
char outputBinFile[] = "a.dat";
char const *arrival_time = "a_";

// compute the mean of an array of samples
float mean(float *samples, int length) {
    float sum = 0.0;

    int i;
    for(i = 0; i < length; i++) {
        sum += samples[i];
    }

    float mean = sum/length;
    return mean;
}

// compute the variance of an array of samples
float variance(float *samples, int length) {
    float smean = mean(samples, length);
    float sum = 0.0;

    int i;
    for(i = 0; i < length; i++) {
        sum += (samples[i]-smean)*(samples[i]-smean);
    }

    float variance = sum/length;
    return variance;
}

// compute the variance of an array of samples assuming 0 mean
float variance_zmean(float *samples, int length) {
    //float smean = mean(samples, length);
    float sum = 0.0;

    int i;
    for(i = 0; i < length; i++) {
        sum += samples[i]*samples[i];
    }

    float variance = sum/length;
    return variance;
}

void separatePS(DataWindow *dataWindow, int separate_win_len, int channel, dtype *pbuffer, dtype *sbuffer) {
   
    int BUFLEN;
    //printf("starting\n");
    int i, j;
    static float a[CHS*CHS];
    int n = CHS;
    int ierr;
    //float *r;
    
    static float w[CHS];
    static float x[CHS*CHS];
    int matz = 1;
    float r;
    float cos_phi;

    int ii;

    static float A[MAXSAMPLINGRATE*(FRONTOFFSET+BACKOFFSET)][CHS];

    static float res[CHS][CHS];

    BUFLEN = dataWindow->size;

    for(i = 0; i < BUFLEN; i++) {
        A[i][0] = dataWindow->buffer0[i];
        A[i][1] = dataWindow->buffer1[i];
        A[i][2] = dataWindow->buffer2[i];
        //printf("%f ", dataWindow->buffer0[i]);
    }

    //printf("0000\n");
    for(ii = 0; ii < BUFLEN; ii++) {
        pbuffer[ii] = A[i][0];
        sbuffer[ii] = A[i][channel];
    }

    //printf("1111\n");
    for(ii = separate_win_len; ii < BUFLEN; ii++) {
        cov(A, ii-separate_win_len, ii, res, BUFLEN);
        for(i = 0; i < CHS; i++) {
            for(j = 0; j < CHS; j++) {
                a[i+j*CHS] = res[i][j];
            }
        }

        ierr = rs ( n, a, w, matz, x );

        if(ierr != 0) {
            printf("!eigenvalue computing error!\n");
        }

        r = 1 - (w[1]+w[0])/(2*w[2]);
        cos_phi = fabs(x[6]);
        pbuffer[ii] = r*cos_phi;
        sbuffer[ii] = r*(1-cos_phi);

    }

    for(i = 0; i < BUFLEN; i++) {
        pbuffer[i] = A[i][0]; //pbuffer[i];
        sbuffer[i] = A[i][channel] * sbuffer[i];
    }

}

// P-phase arrival time picking algorithm
//int picking(dtype *buffer, int size, int pre_len, EventInfo *arrivalTime) {
int picking(dtype *buffer, int size, int pre_len) {

    float sigma_1 = variance(buffer, pre_len);
    int length = size;

    float sigma_2 = 0.0;									//signal variance after the change point
    int k_star = 0;
    float estimator_star = 0.0;

    //printf("buffer0 length: %d\n", length);
    int k;
#ifdef LINUX
    for(k = 0; k < length; k++) {
        buffer[k] = buffer[k];
    }
#else
    for(k = 0; k < length; k++) {
        buffer[k] = buffer[k]*buffer[k];
    }
#endif

    volatile float lnsigma_1 = log(sigma_1);
    float recsigma_1 = 1.0/sigma_1*0.5;

    float sum = 0.0;
    for(k = length-2; k >= 0; k--) {
#ifdef LINUX
        sum += buffer[k+1]*buffer[k+1];
#else
        sum += buffer[k+1];
#endif
        sigma_2 = sum/(length-k);

        volatile float res = log(sigma_2);

        float left_component = 0.5*(lnsigma_1 - res);
        float right_component = 0.5/sigma_2-recsigma_1;
        float estimator = 0.0;

        int snum = length - (k+1);
        estimator = snum*left_component - sum*right_component;

        if(estimator > estimator_star) {
            estimator_star = estimator;
            k_star = k;
        }
    }

    /*
    arrivalTime->ptime = arrivalTime->ptime + (float)(k_star+1)/arrivalTime->sps;

    //test S-wave picking and magnitude
    arrivalTime->stime = arrivalTime->ptime + 2.0;
    arrivalTime->amplitude = 100.0;
    */

    return k_star+1;
}

void getArrivalTime(DataWindow *dataWindow, int separate_win_len, int channel, dtype *pbuffer, dtype *sbuffer, EventInfo *arrivalTime) {

    int BUFLEN;
    BUFLEN = dataWindow->size;

    separatePS(dataWindow, separate_win_len, channel, pbuffer, sbuffer);
    int p_pick, s_pick;
    
    p_pick = picking(pbuffer, BUFLEN, dataWindow->samplingRate);
    s_pick = picking(sbuffer+separate_win_len, BUFLEN-separate_win_len, dataWindow->samplingRate);

    arrivalTime->pTime = arrivalTime->pTime + (double)(p_pick)/dataWindow->samplingRate;
    arrivalTime->sTime = arrivalTime->sTime + (double)(separate_win_len)/dataWindow->samplingRate + (double)(s_pick)/dataWindow->samplingRate;
    arrivalTime->amplitude = maxap(dataWindow->buffer0, BUFLEN);

    // output a test signal
    /*
    int i;
    if(arrivalTime->ptime > 55488.0 && arrivalTime->ptime < 55500.0) {
        FILE *zSig;
        FILE *nSig;
        FILE *eSig;
        FILE *pOut;
        FILE *sOut;

        zSig = fopen("zwave.txt", "w+");
        nSig = fopen("nwave.txt", "w+");
        eSig = fopen("ewave.txt", "w+");
        pOut = fopen("pwave.txt", "w+");
        sOut = fopen("swave.txt", "w+");

        printf("outputing sample data: %d\n", dataWindow->size);
        
        for(i = 0; i < dataWindow->size; i++) {
            fprintf(zSig, "%e\n", dataWindow->buffer0[i]);
            fprintf(nSig, "%e\n", dataWindow->buffer1[i]);
            fprintf(eSig, "%e\n", dataWindow->buffer2[i]);
            fprintf(pOut, "%e\n", pbuffer[i]);
            fprintf(sOut, "%e\n", sbuffer[i]);
        }

        fclose(zSig);
        fclose(nSig);
        fclose(eSig);
        fclose(pOut);
        fclose(sOut);
    }
    */
   
    //end output a test signal

    return;
}

dtype maxap(dtype *buffer, int size) {
    int i;
    dtype max = 0;

    for(i = 0; i < size; i++) {
        if(fabs(buffer[i]) > max)
            max = fabs(buffer[i]);
    }
    return max;
}

BOOL detectEvent(DataWindow *dataWindow, double RATIO_THRESHOLD) {
    static unsigned int rsam_counter = 0;
    static dtype pre_mean = 0.0;
    //static BOOL previous_detected = FALSE;
    static BOOL detected = FALSE;

    static dtype sta_rsams[STA_WIN_SIZE];
    static dtype lta_rsams[LTA_WIN_SIZE];
    static dtype sta_mean;
    static dtype lta_mean;
    static float sta_lta_ratio;

    static unsigned int lta_offset;
    static unsigned int sta_offset;

    dtype rsam_value;
#ifdef LINUX
    float sum = 0.0;
#else
    dtype temp = 0;
    dtype avg = 0;
#endif
   
    rsam_value = 0;
   
    int i;
    //int sps;
    //sps = dataWindow->samplingRate;
    for(i = 0; i < dataWindow->samplingRate; i++) {
#ifdef LINUX
        if(dataWindow->buffer0[i] > 0)
            sum += dataWindow->buffer0[i];
        else
            sum -= dataWindow->buffer0[i];

        if(dataWindow->buffer0[i] > pre_mean)
            rsam_value += (dataWindow->buffer0[i] - pre_mean);
        else
            rsam_value += (pre_mean - dataWindow->buffer0[i]);
#else
        temp = dataWindow->buffer0[i] >> 6;
        if(temp >= pre_mean)
            rsam_value += (temp - pre_mean);
        else
            rsam_value += (pre_mean - temp);

        avg += temp;
#endif
    }

    sta_offset = rsam_counter % STA_WIN_SIZE;
    sta_rsams[sta_offset] = rsam_value;
    lta_offset = rsam_counter % LTA_WIN_SIZE;
    lta_rsams[lta_offset] = rsam_value;

    rsam_counter++;

#ifdef LINUX
    pre_mean = sum/(float)dataWindow->samplingRate;
#else
    pre_mean = avg/dataWindow->samplingRate;
#endif

    if(rsam_counter >= LTA_WIN_SIZE) {
        sta_mean = 0;
        for(i = 0; i < STA_WIN_SIZE; i++) {
            sta_mean += sta_rsams[i];
        }
#ifdef LINUX
        sta_mean = sta_mean/(float)STA_WIN_SIZE;
#else
        sta_mean = sta_mean >> STA_WIN_BIT;
#endif

        lta_mean = 0;
        for(i = 0; i < LTA_WIN_SIZE; i++) {
            lta_mean += lta_rsams[i];
        }
#ifdef LINUX
        lta_mean = lta_mean/(float)LTA_WIN_SIZE;
#else
        lta_mean = lta_mean >> LTA_WIN_BIT;
#endif

        sta_lta_ratio = (float)sta_mean/(float)lta_mean;
    }

    if(sta_lta_ratio > RATIO_THRESHOLD) {
        if(!detected) {
            detected = TRUE;
            return TRUE;
        } else {
            return FALSE;
        }
    } else {
        detected = FALSE;
        return FALSE;
    }
}

/*
void dispEventInfo(EventInfo event_info, int event_counter){
    printf("\nevent %d:\n  Long = %f;\n  Lat = %f;\n  Alt = %f;\n  PTime = %f;\n  STime = %f;\n  Amp = %f\n", event_counter, event_info.x, event_info.y,event_info.z,event_info.ptime,event_info.stime,event_info.amplitude);

//    printf("\nevent %d: 00000: 0000: %f\n", event_counter, event_info.ptime);
}
*/

int dataProcess(BufferInput *bufInput, EventInfo *event_info, int id, double *streamStartTime, dtype *streamBuffer, double threshold) {
     
    static dtype pbuffer[MAXSAMPLINGRATE*(FRONTOFFSET+BACKOFFSET)], sbuffer[MAXSAMPLINGRATE*(FRONTOFFSET+BACKOFFSET)];
    static dtype *buffer[3];

    // test for picking
    //static double lastPicking = 0.0;
    //static int wPicking = 0;

    // count total events since the first process;
    static int totalEventCounter = 0;

    static BOOL newstart = TRUE;
    int extraoffset;
    int i, k;
    double starttime, timingoffset;
    int samplingrate;
    int numseconds;
    
    samplingrate = bufInput->samplingRate;
    numseconds = bufInput->size/bufInput->samplingRate;
    timingoffset = bufInput->timingOffset;

    if(newstart == TRUE) {
        extraoffset = BACKOFFSET;
        starttime = timingoffset;

        for(k = 0; k < CHS; k++) {
            buffer[k] = (dtype*)malloc((numseconds+FRONTOFFSET+BACKOFFSET)*samplingrate*sizeof(dtype));
        }
        newstart = FALSE;
    } else {
        //printf("false timing offset\n");
        extraoffset = 0; 
        starttime = timingoffset-BACKOFFSET;
    }

    //printf("extraoffset: %d\n", extraoffset);

    
    for(k = 0; k < CHS; k++) {
        memcpy(buffer[k]+(FRONTOFFSET+BACKOFFSET)*samplingrate, bufInput->buffer+k*samplingrate*numseconds, bufInput->size*sizeof(dtype));
    }

    numseconds += FRONTOFFSET+BACKOFFSET; //give total buffer size including front and back offset
    
    DataWindow data_win;
    
	int frontoffset, backoffset;
    
    BOOL event_detected;
	event_detected = FALSE;

	frontoffset = FRONTOFFSET + extraoffset, backoffset = BACKOFFSET;

    int event_counter;
	event_counter = 0;
    
    data_win.samplingRate = samplingrate;

    for(i = frontoffset; i < numseconds-backoffset; i++) {
        data_win.buffer0 = buffer[0] + i*samplingrate;
        data_win.buffer1 = buffer[1] + i*samplingrate;
        data_win.buffer2 = buffer[2] + i*samplingrate;
            
        event_detected = detectEvent(&data_win, threshold);
        
        if(event_detected) {
            data_win.buffer0 = data_win.buffer0 - FRONTOFFSET*samplingrate;
            data_win.buffer1 = data_win.buffer1 - FRONTOFFSET*samplingrate;
            data_win.buffer2 = data_win.buffer2 - FRONTOFFSET*samplingrate;
            data_win.size = (FRONTOFFSET+BACKOFFSET)*samplingrate;
            
            // fill up the event_info
            event_info[event_counter].stationID = bufInput->stationID;
            event_info[event_counter].samplingRate = samplingrate;
            event_info[event_counter].x = bufInput->x;
            event_info[event_counter].y = bufInput->y;
            event_info[event_counter].z = bufInput->z;
            event_info[event_counter].pTime = starttime + (i - frontoffset) - FRONTOFFSET;
            event_info[event_counter].sTime = starttime + (i - frontoffset) - FRONTOFFSET;
                
            streamStartTime[event_counter] = event_info[event_counter].pTime;
            getArrivalTime(&data_win, PSSEPARATEWIN*samplingrate, 1, pbuffer, sbuffer, &event_info[event_counter]);

            // memcpy(streamBuffer+(FRONTOFFSET+BACKOFFSET)*samplingrate*event_counter, pbuffer, samplingrate*(FRONTOFFSET+BACKOFFSET)*sizeof(dtype));
            // printf("event %d: %d: 0000: %f: %f\n", event_counter, i-FRONTOFFSET, event_info[event_counter].ptime, event_info[event_counter].stime);

            // downsampling stream data for sending out
            size_t jOffset = (FRONTOFFSET+BACKOFFSET)*samplingrate*event_counter/DOWN_SAMPLING_FAC;
            size_t jLength = (FRONTOFFSET+BACKOFFSET)*samplingrate/DOWN_SAMPLING_FAC;
            //printf("+++++++++ prepare sending buffer %zu, %zu\n", jOffset, (jOffset+jLength));
            
            for(size_t j = 0; j < jLength; j++) {
                streamBuffer[j+jOffset] = pbuffer[j*DOWN_SAMPLING_FAC];
            }
            
            event_counter++;
            i += BACKOFFSET;
        }
    }

    for(k = 0; k < 3; k++) { // prepare space for next file
        memcpy(buffer[k], buffer[k] + (numseconds - FRONTOFFSET - BACKOFFSET)*samplingrate, (FRONTOFFSET+BACKOFFSET)*samplingrate*sizeof(dtype));
    }

    if(event_counter > 0) {
        for(k = 0; k < event_counter; k++) {
            totalEventCounter++;

            //if(event_info[k].pTime < lastPicking) {
            //    wPicking++;
            //    printf("%d >>>>>>>>>>>>>>>>>>>>\n", wPicking);
            //}

            //printf("    %d, %lf\n", totalEventCounter, event_info[k].pTime);
            //lastPicking = event_info[k].pTime;
            
            WriteEventTxt(outputTxtFile, &event_info[k], totalEventCounter);
            WriteEventBin(outputBinFile, &event_info[k], totalEventCounter, id);
        }
    }

	return event_counter;
}

// convert day of year to date (year/month/day)
Date DayOfYearToDate(int year, int dayOfYear) {
    Date date;

    static const int MonthLen[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    int leap;
    int month;
    int dayOfMonth;

    leap = (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);
    dayOfMonth = dayOfYear;

    for (month = 0; month < 12; month++) {
        int mLen = MonthLen[month];
        if (leap && month == 1)
            mLen++;
        if (dayOfMonth <= mLen)
            break;

        dayOfMonth -= mLen;
    }

    date.year = year;
    date.month = month;
    date.dayOfMonth = dayOfMonth;

    return date;
}

// convert date/time to epoch time
double ToEpochTime(int year, int month, int day, int hour, int minute, int second, int milliSecond, int dstOn) {

    double epochDouble;

    struct tm t;
    time_t t_of_day;
    
    t.tm_year = year - 1900;
    t.tm_mon = month - 1;     // Month, 0 - jan
    t.tm_mday = day;          // Day of the month
    t.tm_hour = hour;     
    t.tm_min = minute;
    t.tm_sec = second;
    t.tm_isdst = -1;          // Is DST on? 1 = yes, 0 = no, -1 = unknown
    t_of_day = mktime(&t);

    // printf("%d-%d-%d %d:%d:%d:%d\n", year, month, day, hour, minute, second, milliSecond);
    // printf("-------- %ld\n", t_of_day); 

    // printf("The current local time is: %s\n", ctime(&t_of_day));

    epochDouble = (double)t_of_day;
    epochDouble += (double)milliSecond/1000.0;

    return epochDouble;
}

void itoa(char *buf, int base, int d) {
    char *p = buf;
    char *p1, *p2;
    unsigned long ud = d;
    int divisor = 10;

    /* If %d is specified and D is minus, put `-' in the head.  */
    if (base == 'd' && d < 0) {
        *p++ = '-';
        buf++;
        ud = -d;
    } else if (base == 'x') {
        divisor = 16;
    }

    /* Divide UD by DIVISOR until UD == 0.  */
    do {
        int remainder = ud % divisor;

        *p++ = (remainder < 10) ? remainder + '0' : remainder + 'a' - 10;
    } while (ud /= divisor);

    /* Terminate BUF.  */
    *p = 0;

    /* Reverse BUF.  */
    p1 = buf;
    p2 = p - 1;
    while (p1 < p2) {
        char tmp = *p1;
        *p1 = *p2;
        *p2 = tmp;
        p1++;
        p2--;
    }
}

// write arrival time information to file
void WriteEventTxt(char *arrivalTimeFile, EventInfo *event, int eventCounter) {
    FILE *fptr;
    struct tm *pTime;
    struct tm *sTime;

    long pTimeInSec;
    long sTimeInSec;

    time_t pTime_t;
    time_t sTime_t;

    double pMilliSec;
    double sMilliSec;

    pTimeInSec = (long) event->pTime;
    sTimeInSec = (long) event->sTime;

    pMilliSec = event->pTime - (double)pTimeInSec;
    sMilliSec = event->sTime - (double)sTimeInSec;

    pTime_t = (time_t)pTimeInSec;
    sTime_t = (time_t)sTimeInSec;

    fptr = fopen(arrivalTimeFile, "a+");
    if(fptr != NULL) {

        event->eventSeqNum = eventCounter;

        pTime = localtime(&pTime_t);

        fprintf(fptr, "%04d,%05d,%4d,%lf,%lf,%lf,%02d/%02d/%04d %02d:%02d:%06.3lf,", 
                event->eventSeqNum, event->stationID, event->samplingRate, event->x, event->y, event->z, 
                pTime->tm_mon+1, pTime->tm_mday, pTime->tm_year+1900, pTime->tm_hour, 
                pTime->tm_min, (double)pTime->tm_sec+pMilliSec);
        
        sTime = localtime(&sTime_t);
        
        fprintf(fptr, "%02d/%02d/%04d %02d:%02d:%06.3lf,%.12lf\n", 
                sTime->tm_mon+1, sTime->tm_mday, sTime->tm_year+1900, sTime->tm_hour, 
                sTime->tm_min, (double)sTime->tm_sec+sMilliSec, event->amplitude);
    }

    fclose(fptr);
}

void WriteEventBin(char *arrivalTimeFile, EventInfo *event, int eventCounter, int id) {

    FILE *fptr;
    char id_a[10],arrival_timefile[80];
    itoa(id_a, 10, id);
    strcpy(arrival_timefile, arrival_time);
	strcat(arrival_timefile, id_a);
	strcat(arrival_timefile, ".dat");
    fptr = fopen(arrival_timefile, "ab+");
    if(fptr != NULL) {
        event->eventSeqNum = eventCounter;
        fwrite(event, sizeof(EventInfo), 1, fptr);
    }

    fclose(fptr);
    
}

// convert decimal degrees to seconds
int DegToSec(double decimalDegree) {
    int sign;
    int hour, minute, second;
    int totalSeconds;

    sign = (decimalDegree > 0) ? 1 : -1;
    // printf("sign: %d\n", sign);
    decimalDegree = decimalDegree * ((double)sign);

    hour = (int)decimalDegree;
    minute = (int)((decimalDegree - (double)hour)*60.0);
    second = (int)((decimalDegree - (double)hour - (double)minute/60.0)*3600.0 + 0.5);

    // printf("%d, %d, %d, %d\n", sign, hour, minute, second);

    totalSeconds = sign * (hour*3600 + minute*60 + second);

    return totalSeconds;
}

// convert seconds to decimal degrees
double SecToDeg(int totalSeconds) {
    double sign;
    int hour, minute, second;
    double decimalDegree;

    sign = (totalSeconds >=0 ) ? 1.0 : -1.0;

    totalSeconds = totalSeconds * ((int)sign);
    
    hour = totalSeconds / 3600;
    minute = (totalSeconds - hour*3600) / 60;
    second = totalSeconds - hour*3600 - minute*60;

    // printf(".... %d %d %d\n", hour, minute, second);
    decimalDegree = sign * ((double)hour + (double)minute/60.0 + (double)second/3600.0);

    return decimalDegree;
}
