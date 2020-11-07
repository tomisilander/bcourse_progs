#include <signal.h>
#include <unistd.h>

int rep_flag_set;
int end_flag_set;
int alrm_flag_set;

#define SIGHANDLER(SIGNAME, FLAGNAME)\
void SIGNAME ## _handler(int signum) {\
  signal(signum, SIGNAME ## _handler);\
  signum = 0;\
  FLAGNAME ## _set = 1;\
}

SIGHANDLER(usr1, rep_flag);
SIGHANDLER(usr2, end_flag);
SIGHANDLER(alrm, alrm_flag);

void set_signal_handlers(){
  rep_flag_set  = 0;
  end_flag_set  = 0;
  alrm_flag_set = 0;
  signal(SIGUSR1, usr1_handler);
  signal(SIGUSR2, usr2_handler);
  signal(SIGALRM, alrm_handler);
}

void setalarm(int s) {
  if(s > 0) alarm(s);
}
