#define _group_represent(v,u)\
 (v).c[0] = -(-cimag((u).c[0])*cimag((u).c[3])+cimag((u).c[1])*cimag((u).c[2])-creal((u).c[0])*creal((u).c[3])+creal((u).c[1])*creal((u).c[2]));\
 (v).c[1] = -(-cimag((u).c[0])*creal((u).c[3])-cimag((u).c[1])*creal((u).c[2])+cimag((u).c[2])*creal((u).c[1])+cimag((u).c[3])*creal((u).c[0]));\
 (v).c[2] = +cimag((u).c[0])*creal((u).c[2])-cimag((u).c[1])*creal((u).c[3])-cimag((u).c[2])*creal((u).c[0])+cimag((u).c[3])*creal((u).c[1]);\
 (v).c[3] = -cimag((u).c[0])*creal((u).c[3])+cimag((u).c[1])*creal((u).c[2])-cimag((u).c[2])*creal((u).c[1])+cimag((u).c[3])*creal((u).c[0]);\
 (v).c[4] = +cimag((u).c[0])*cimag((u).c[3])+cimag((u).c[1])*cimag((u).c[2])+creal((u).c[0])*creal((u).c[3])+creal((u).c[1])*creal((u).c[2]);\
 (v).c[5] = -(-cimag((u).c[0])*cimag((u).c[2])+cimag((u).c[1])*cimag((u).c[3])-creal((u).c[0])*creal((u).c[2])+creal((u).c[1])*creal((u).c[3]));\
 (v).c[6] = -(+cimag((u).c[0])*creal((u).c[1])-cimag((u).c[1])*creal((u).c[0])-cimag((u).c[2])*creal((u).c[3])+cimag((u).c[3])*creal((u).c[2]));\
 (v).c[7] = -(-cimag((u).c[0])*cimag((u).c[1])+cimag((u).c[2])*cimag((u).c[3])-creal((u).c[0])*creal((u).c[1])+creal((u).c[2])*creal((u).c[3]));\
 (v).c[8] = +4.999999999999999e-01*(+cimag((u).c[0])*cimag((u).c[0])-cimag((u).c[1])*cimag((u).c[1])-cimag((u).c[2])*cimag((u).c[2])+cimag((u).c[3])*cimag((u).c[3])+creal((u).c[0])*creal((u).c[0])-creal((u).c[1])*creal((u).c[1])-creal((u).c[2])*creal((u).c[2])+creal((u).c[3])*creal((u).c[3]));\


