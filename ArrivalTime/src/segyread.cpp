// 
// apis to read segy files header and data
//
// Author: Lei Shi
//

#include <stdio.h>
#include <string.h>
#include "segy.h"    /* SEGY header structure */
#include "ArrivalTime.h"

// get SEGY file trace header information
SEGY_IO_STATUS segy_head_read(const char *filename, SEGYTRACEHEAD *header) {
    
    FILE *input;
    input = fopen(filename, "r");

    if(input == NULL) {
        return FILE_OPEN_FAILED;
    } 

    size_t headerSize = fread(header, sizeof(SEGYTRACEHEAD), 1, input);
   

    if(headerSize == 0) {
        fclose(input);
        return FILE_OPEN_FAILED;
    }
    
    fclose(input);

    return SUCCESS;
}

// get SEGY file data and put it in a buffer
SEGY_IO_STATUS segy_read(const char *filename, SEGYTRACEHEAD *header, dtype *buffer, int size) {

    //printf("++++++++ %s\n\n", filename);
    FILE *input;
    input = fopen(filename, "r");

    if(input == NULL) {
        fclose(input);
        return FILE_OPEN_FAILED;
    } 

    size_t headerSize = fread(header, sizeof(SEGYTRACEHEAD), 1, input);
    if(headerSize == 0) {
        fclose(input);
        return FILE_OPEN_FAILED;
    }

    //printf("    getting header....\n");

    int ele_size = 0; // element size
    int element;

    //int temp;
    //unsigned int sample;

    if(header->data_form == 0) {
        ele_size = 2;
    } else if(header->data_form == 1) {
        ele_size = 4;
    }

    //printf("    getting form: %d\n", header->data_form);

    int buf_pt = 0;

    //unsigned int max = 0;
    //unsigned int min = 0;

    size_t eleSize;
    while(!feof(input)) {

        eleSize = fread(&element, ele_size, 1, input);
        
        if(buf_pt < size && eleSize == 1) {
            buffer[buf_pt] = (dtype)element * header->scale_fac;
        }

        buf_pt++;
    }


//////////////////////////////////////
 /*   char* str1; char* str2;
str1 = "data.dat";
//(char *)malloc(strlen(filename)+strlen(str1));

     str2 = strcat((char *)filename,str1);

//printf("++++++++ %s\n", str2);
FILE *fp;
fp=fopen(str2,"w");
fwrite(buffer,ele_size,size,fp);
fclose(fp);
*/
    //printf("++++++++ %u, %u\n", min, max);

    fclose(input);

    return SUCCESS;
}

