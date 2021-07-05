
#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H
typedef enum
{
    STORE,
    DONTSTORE
} storage_switch;

typedef struct _data_storage
{
    int n;
    int *ni;
    double *data;
} data_storage;

typedef struct _data_storage_array
{
    int n;
    data_storage *element;

} data_storage_array;





#endif