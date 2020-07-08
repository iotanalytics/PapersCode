#ifndef SEGY_READ
#define SEGY_READ

#include "segy.h"

SEGY_IO_STATUS segy_head_read(const char *filename, SEGYTRACEHEAD *header);

SEGY_IO_STATUS segy_read(const char *filename, SEGYTRACEHEAD *header, dtype *buffer, int size);

#endif
