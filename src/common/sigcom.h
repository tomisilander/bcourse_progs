#ifndef SIGCOM_H
#define SIGCOM_H

extern int rep_flag_set;
extern int end_flag_set;
extern int alrm_flag_set;

extern void set_signal_handlers();
extern void setalarm(int s);

#endif
