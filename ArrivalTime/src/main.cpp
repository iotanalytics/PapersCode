#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ArrivalTime.h"
#include "segyread.h"

extern "C" {
#include "dbh.h" 
}
//#include "socket.h"

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <sys/stat.h>
#include <sys/signal.h>
#include <dirent.h>

#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <pthread.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <termios.h>
#include <unistd.h>
#include <time.h>
#include <cstring>
#include <math.h>

using namespace std;
//#### test locally without touching
int xbee_fd;
unsigned char stream_control;
pthread_t listener_thread,server_thread,reader_thread;
unsigned char dest_address[8];

double threshold = 2.0;
int filter = 0;
int low_freq = 2;
int high_freq = 10;

int list_files(string curdir, vector< vector<string> > &file_lists) {
    
    vector<string> file_list;
    DIR *d;
    struct dirent *dir;

    for(int i = 0; i < 3; i++) {
        file_lists.push_back(file_list);
    }

    cout << "open dir " << curdir << endl;
    d = opendir(curdir.c_str());

    size_t found;
    string extz("HZ");
    string extn("HN");
    string exte("HE");
    //if (found != string::npos) {
    //    cout << "first 'needle' found at: " << found << '\n';
    //}

    if (d) {
        while ((dir = readdir(d)) != NULL) {
            string file_name(dir->d_name);
            //cout << "r file name: " << file_name << endl;
            found = file_name.find(extz);
            if(found != string::npos) {
                string full_path = curdir + "/" + file_name;
                file_lists[0].push_back(full_path);
            }
            found = file_name.find(extn);
            if(found != string::npos) {
                string full_path = curdir + "/" + file_name;
                file_lists[1].push_back(full_path);
            }
            found = file_name.find(exte);
            if(found != string::npos) {
                string full_path = curdir + "/" + file_name;
                file_lists[2].push_back(full_path);
            }
        }

        closedir(d);
    }

    sort( file_lists[0].begin(), file_lists[0].end() );
    sort( file_lists[1].begin(), file_lists[1].end() );
    sort( file_lists[2].begin(), file_lists[2].end() );

    /*
    for(int i = 0; i < 3; i++) {
        cout << "file list " << i << endl;
        for(size_t j = 0; j < file_lists[i].size(); j++) {
            cout << "    " << file_lists[i][j] << endl;
        }
    }*/

    int file_list_length = file_lists[0].size();
    return file_list_length;
}

void ConfigParse(const char *configFile, map<string, string> &configs) {

    ifstream inFile(configFile, ifstream::in);

    string line;
    while(getline(inFile, line)) {
        istringstream isLine(line);
        string key;

        if(getline(isLine, key, '=')) {
            string value;
            if(getline(isLine, value)) {
                configs[key] = value;

                cout << key << " - - - - " << value << endl;
            }
        }
    }

}

/*
//test GUI monitoring
void SendArrivalToGUI(EventInfo *event) {
   
    unsigned char src_add[3];
    struct eventMsgGUI{
        unsigned char src_add[3];		
		unsigned int length;
		unsigned int type;		
		EventInfo event;	
	} msg;

    unsigned char arrtimeBuffer[58];

    arrtimeBuffer[0]=DP_TYPE;
    memcpy(arrtimeBuffer+1,src_add,3);
    msg.type=DP_TYPE;	
    
	msg.event=*event;
	msg.length=DP_PLD_LEN;
    memcpy(arrtimeBuffer+4,&msg.length,2);	
    memcpy(arrtimeBuffer+6,&event->stationID,2);
    memcpy(arrtimeBuffer+8,&event->samplingRate,2);
    memcpy(arrtimeBuffer+10,&event->x,8);
    memcpy(arrtimeBuffer+18,&event->y,8);
    memcpy(arrtimeBuffer+26,&event->z,8);
    memcpy(arrtimeBuffer+34,&event->pTime,8);
    memcpy(arrtimeBuffer+42,&event->sTime,8);
    memcpy(arrtimeBuffer+50,&event->amplitude,8);  
		
    int val=sendToClient(arrtimeBuffer,58);
	printf("********Arrtime information Broadcasted to %d clients*********\n",val);
}

void SendStreamToGUI() {

}
*/

int main(int argc, char *argv[]) {
  
    // printf("number of arguments: %d\n", argc);

    if(argc != 4) {
        printf("Usage: %s [input_seismic_files_folder] [station_id] [config_file]\n", argv[0]); 
        exit(1);
    }

    // read config file
    map<string, string> configs;
    map<string, string>::iterator ssit;

    cout << "parsing" << endl;
    ConfigParse(argv[3], configs);

    for(ssit = configs.begin(); ssit != configs.end(); ++ssit) {
        cout << ssit->first << ": " << ssit->second << endl;
    }
    cout << "parse done" << endl;

    stringstream configValue;

    configValue.clear();
    configValue.str(configs["threshold"]);
    configValue >> threshold;
    cout << "threshold: " << threshold << endl;
   
    configValue.clear();
    configValue.str(configs["filter"]);
    configValue >> filter;
    cout << "filter: " << filter << endl;
 
    configValue.clear();
    configValue.str(configs["low_freq"]);
    configValue >> low_freq;
    cout << "low_freq: " << low_freq << endl;

    configValue.clear();
    configValue.str(configs["high_freq"]);
    configValue >> high_freq;
    cout << "high_freq: " << high_freq << endl;

    /*
    if(!pthread_create(&listener_thread, NULL, &tcp_socket_client, NULL)) {
	printf("Listener/tcp client spawned successfully\n");
	    arg.pid=listener_thread;
	    fcntl(xbee_fd,15,&arg);
    } else {
	    printf("Error creating server thread\n");
	    return 1;
    }*/
    
    /*
    if(!pthread_create(&server_thread, NULL, &tcp_socket_server, NULL)) {
	    printf("TCP server spawned successfully\n");
    } else {
        printf("Error creating server thread\n");
	    return 1;
    }*/

    //while(1) {}
    //sleep(20);

    SEGYTRACEHEAD segyheader;        // segy header of each file
    //dtype *buffer[3];             // data buffer0

    //FILE *station_file;
    //FILE *arrival_time_file;
    
    int samplingrate, k, station_id;
    int num_samples;
    int year;
    int month;
    int day;
    int dayOfYear;
    int hour;
    int min;
    int sec;
    int millisec;
    double timingoffset;
    int numseconds;
    string ftp_map_directory = string(argv[1]);

    // #### test locally without touching
    vector< vector<string> > file_lists;
    int n_files = list_files(ftp_map_directory, file_lists);
    printf("total number of files: %d\n", n_files);

    // arrival_time_file = fopen(argv[2], "ab+");
    station_id = atoi(argv[2]);
    printf("\n\n\n\nwe have station %d,\n", station_id);

    /*
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        cout << ' ' << *it;
        cout << '\n';
    }*/
    // #### end
    
    //vector<EventInfo> eventinfos;
    //vector<SEGYTRACEHEAD> segyheaders;

    // load segy file header
    SEGY_IO_STATUS status;
    int event_counter=0, total_counter=0;

    BufferInput bufInput;
    EventInfo event_info[MAXNUMEVENTS];
    double streamStartTime[MAXNUMEVENTS];
    dtype *streamBuffer;
    streamBuffer = (dtype*)malloc((FRONTOFFSET+BACKOFFSET)*MAXSAMPLINGRATE*sizeof(dtype)*MAXNUMEVENTS/DOWN_SAMPLING_FAC);

    for(size_t i = 0; i < file_lists[0].size(); i++) { 
        // Load segy header information
        for(k=0; k<CHS; k++) {
            status = segy_head_read(file_lists[k][i].c_str(), &segyheader); 
            if(status == FILE_OPEN_FAILED) {
                printf("File open file %s failed!\n", file_lists[k][i].c_str());
                exit(1);
	    	} else {
                if(0 == k) {
                    //segyheaders.push_back(segyheader);
                }
            }
        }
     
        year = segyheader.year;
        dayOfYear = segyheader.day;
	    hour = segyheader.hour;
	    min = segyheader.minute;
	    sec = segyheader.second;
	    millisec = segyheader.m_secs;
	    samplingrate = (int)(1000.0/((double)segyheader.deltaSample/1000.0));
	    //timingoffset = (double)hour*3600.0 + (double)min*60.0 + (double)sec + (double)millisec/1000;
	    num_samples = segyheader.num_samps;
	    numseconds = num_samples/samplingrate;
	   
        //printf("        %d-%d\n", year, dayOfYear);
        Date date = DayOfYearToDate(year, dayOfYear);
        month = date.month;
        day = date.dayOfMonth;
        //printf("        %d-%d-%d\n", year, month, day);

        timingoffset = ToEpochTime(year, month, day, hour, min, sec, millisec, -1);

        //printf("timing offset    : %lf\n", timingoffset);
        //printf("number of samples: %d\n", num_samples);
        //printf("sample length    : %d\n", segyheader.sampleLength);
        //printf("sampling rate    : %d\n", samplingrate);
        //printf("number of seconds: %d\n", numseconds);

        if (i == 0) {
            bufInput.channelNum = CHS;
            bufInput.buffer = (dtype*)malloc(numseconds*samplingrate*sizeof(dtype)*bufInput.channelNum);
        }

        //printf(">>>>>>>> malloc done!\n");
        
        for(k = 0; k < CHS; k++){
            status = segy_read(file_lists[k][i].c_str(), &segyheader, bufInput.buffer+k*numseconds*samplingrate, num_samples);
        }
        
        bufInput.samplingRate = samplingrate;
        bufInput.stationID = station_id;
        bufInput.size = numseconds*samplingrate;
        bufInput.timingOffset = timingoffset;

        //printf("long: %d, lat: %d, z: %d\n", segyheader.recLongOrX, segyheader.recLatOrY, segyheader.recElevation);
        //printf("Unit: %d, Scale: %d\n", segyheader.coordUnits, segyheader.coordScale);

        bufInput.x = SecToDeg(segyheader.recLongOrX);
        bufInput.y = SecToDeg(segyheader.recLatOrY);
        bufInput.z = (double)segyheader.recElevation;

        //printf("x: %lf, y: %lf\n, z: %lf\n", bufInput.x, bufInput.y, bufInput.z);

        //apply bandpass filter
        
       if(1 == filter) {
           xapiir_rmean(bufInput.buffer, numseconds*samplingrate*bufInput.channelNum, low_freq, high_freq, (double)segyheader.deltaSample/1000000.0); // the last parameter is to convert macrosecond to second
        }
        
        event_counter = dataProcess(&bufInput, event_info, station_id, streamStartTime, streamBuffer, threshold);

        // test stream buffer
        // FILE * streamFile;
        if(event_counter != 0) {

            printf("%d events:\n", event_counter);
            for(int events = 0; events < event_counter; events++) {
                total_counter++;

                /*
                stringstream ss;
                string eFile;
                ss << total_counter;
                eFile = ss.str();
                eFile += ".event";

                streamFile = fopen(eFile.c_str(), "w");
                for(int aa = events*samplingrate*(FRONTOFFSET+BACKOFFSET)/DOWN_SAMPLING_FAC; 
                        aa < (events+1)*samplingrate*(FRONTOFFSET+BACKOFFSET)/DOWN_SAMPLING_FAC; 
                        aa++) {
                    fprintf(streamFile, "%.16f\n", (float)streamBuffer[aa]);
                }

                fclose(streamFile);
                */

                //printf("++++++++ %s\n", bufInput);
                printf("%d, %lf, %lf\n", total_counter, event_info[events].pTime, event_info[events].sTime);

            }

            //total_counter += event_counter;
        }

    } //end  for(int i = 0; i < file_lists[0].size(); i++)

    printf("total event picked: %d\n", total_counter);

    /*
    double a = 0.000000873281934;
    double b = 3927138.3108783748;

    printf("== %.16lf, %.16lf\n", a, b);

    float c = (float)a;
    float d = (float)b;

    printf("== %.16f, %.16f\n", c, d);
    */

    free(bufInput.buffer);
    free(streamBuffer);
    
    return 0;
}


