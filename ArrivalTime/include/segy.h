#ifndef SEGY_H
#define SEGY_H

#define SEGY_EBCDIC_HEADER_SIZE 3200
#define SEGY_BINARY_HEADER_SIZE 400
#define SEGY_TRACE_HEADER_SIZE 240

typedef struct SegyEBCDICHeader {
    char StringHeader[SEGY_EBCDIC_HEADER_SIZE];
} SEGYEBCDICHEADER;

typedef union SegyBinaryHeader {
    char cdata[SEGY_BINARY_HEADER_SIZE];
    struct S400B { //used ones for tape header, follows 3200 header
        int   iJobIdNumber;               //1-4
        int   iLineNumber;                //5-8
        int   iReelNumber;                //9-12
        short siNumDataTracesPerEnsemble; //13-14
        short siNumAuxTracesPerEnsemble;  //15-16
        short siSampleIntInMicro;         //17-18
        short siSampleIntOrig;            //19-20
        short siNumSamplesPerTrace;       //21-22
        short siNumSamplesOrig;           //23-24
        short siDataSampleFormatCode;     //25-26
        short siEnsembleFold;             //27-28
        short siTraceSortingCode;         //29-30
        short siVerticalSumCode;          //31-32
        short siSweepStartFreq;           //33-34
        short siSweepEndFreq;             //35-36
        short siSweepLenMills;            //37-38
        short siSweepTypeCode;            //39-40
        short siSweepChanNumber;          //41-42
        short siSweepTaperLenStart;       //43-44
        short siSweepTaperLenEnd;         //45-46
        short siSweepTaperTypeCode;       //47-48
        short siCorrDataTracesCode;       //49-50
        short siBinaryGainRecovered;      //51-52
        short siAmplitudeRecoveryMethod;  //53-54
        short siMeasurementSysCode;       //55-56
        short siImpulseSigPolarity;       //57-58
        short siVibratoryPolarityCode;    //59-60
        short siChunk;                    //61-62
        int   uiNumSamplesLong;           //63-66
        char  cUndefined1[234];           //67-300
        short siSEGYRevision;             //301-302
        short siFixedLengthTraceFlag;     //303-304
        short siExtendedHeaderCount;      //305-306
  } HEAD400B;
} SETYBINARYHEADER;

typedef struct SegyTraceHead {   /*  Offset Description  */
    int  lineSeq;           /*   0 Sequence numbers within line */
    int  reelSeq;           /*   4 Sequence numbers within reel */
    int  event_number;      /*   8 Original field record number or trigger number */
    int  channel_number;    /*  12 Trace channel number within the original field record */
    int  energySourcePt;    /*  16 X */
    int  cdpEns;            /*  20 X */
    int  traceInEnsemble;   /*  24 X */
    short traceID;           /*  28 Trace identification code: seismic data = 1 */
    short vertSum;           /*  30 X */
    short horSum;            /*  32 X */
    short dataUse;           /*  34 X */
    int  sourceToRecDist;   /*  36 X */
    int  recElevation;      /*  40 X */
    int  sourceSurfaceElevation; /*  44 X */
    int  sourceDepth;       /*  48 X */
    int  datumElevRec;      /*  52 X */
    int  datumElevSource;   /*  56 X */
    int  sourceWaterDepth;  /*  60 X */
    int  recWaterDepth;     /*  64 X */
    short elevationScale;    /*  68 Elevation Scaler: scale = 1 */
    short coordScale;        /*  70 Coordinate Scaler: scale = 1 */
    int  sourceLongOrX;     /*  72 X */
    int  sourceLatOrY;      /*  76 X */
    int  recLongOrX;        /*  80 X */
    int  recLatOrY;         /*  84 X */
    short coordUnits;        /*  88 Coordinate Units:  = 2 (Lat/Long) */
    short weatheringVelocity;/*  90 X */
    short subWeatheringVelocity; /*  92 X */
    short sourceUpholeTime;  /*  94 X */
    short recUpholeTime;     /*  96 X */
    short sourceStaticCor;   /*  98 X */
    short recStaticCor;      /* 100 X */
    short totalStatic;       /* 102 Total Static in MILLISECS added to Trace Start Time (lower 2 bytes)*/
    short lagTimeA;          /* 104 X */
    short lagTimeB;          /* 106 X */
    short delay;             /* 108 X */
    short muteStart;         /* 110 X */
    short muteEnd;           /* 112 X */
    short sampleLength;      /* 114 Number of samples in this trace (unless == 32767) */
    short deltaSample;       /* 116 Sampling interval in MICROSECONDS (unless == 1) */
    short gainType;          /* 118 Gain Type: 1 = Fixed Gain */
    short gainConst;         /* 120 Gain of amplifier */
    short initialGain;       /* 122 X */
    short correlated;        /* 124 X */
    short sweepStart;        /* 126 X */
    short sweepEnd;          /* 128 X */
    short sweepLength;       /* 130 X */
    short sweepType;         /* 132 X */
    short sweepTaperAtStart; /* 134 X */
    short sweepTaperAtEnd;   /* 136 X */
    short taperType;         /* 138 X */
    short aliasFreq;         /* 140 X */
    short aliasSlope;        /* 142 X */
    short notchFreq;         /* 144 X */
    short notchSlope;        /* 146 X */
    short lowCutFreq;        /* 148 X */
    short hiCutFreq;         /* 150 X */
    short lowCutSlope;       /* 152 X */
    short hiCutSlope;        /* 154 X */
    short year;              /* 156 year of Start of trace */
    short day;               /* 158 day of year at Start of trace */
    short hour;              /* 160 hour of day at Start of trace */
    short minute;            /* 162 minute of hour at Start of trace */
    short second;            /* 164 second of minute at Start of trace */
    short timeBasisCode;     /* 166 Time basis code: 2 = GMT */
    short traceWeightingFactor; /* 168 X */
    short phoneRollPos1;     /* 170 X */
    short phoneFirstTrace;   /* 172 X */
    short phoneLastTrace;    /* 174 X */
    short gapSize;           /* 176 X */
    short taperOvertravel;   /* 178 X */
    char  station_name[6];   /* 180 Station Name code (5 chars + \0) */
    char  sensor_serial[8];  /* 186 Sensor Serial code (7 chars + \0) */
    char  channel_name[4];   /* 194 Channel Name code (3 chars + \0) */
    short totalStaticHi;     /* 198 Total Static in MILLISECS added to Trace Start Time (high 2 bytes)*/
    int  samp_rate;         /* 200 Sample interval in MICROSECS as a 32 bit integer */
    short data_form;         /* 204 Data Format flag: 0=16 bit, 1=32 bit integer */
    short m_secs;            /* 206 MILLISECONDS of seconds of Start of trace */
    short trigyear;          /* 208 year of Trigger time */
    short trigday;           /* 210 day of year at Trigger time */
    short trighour;          /* 212 hour of day at Trigger time */
    short trigminute;        /* 214 minute of hour at Trigger time */
    short trigsecond;        /* 216 second of minute at Trigger time */
    short trigmills;         /* 218 MILLISECONDS of seconds of Trigger time */
    float scale_fac;         /* 220 Scale Factor (IEEE 32 bit float) */
    short inst_no;           /* 224 Instrument Serial Number */
    short not_to_be_used;    /* 226 X */
    int  num_samps;         /* 228 Number of Samples as a 32 bit integer
                              * (when sampleLength == 32767) */
    int  max;               /* 232 Maximum value in Counts */
    int  min;               /* 236 Minimum value in Counts */
} SEGYTRACEHEAD;                 /* end of segy trace header */

typedef enum { SUCCESS, FILE_OPEN_FAILED } SEGY_IO_STATUS;

#endif                      /* SEGY_H        */
