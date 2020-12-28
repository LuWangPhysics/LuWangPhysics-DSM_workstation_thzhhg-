#include <time.h>
struct mytime
{       int t1;
	int t2;
	mytime(){}
	~mytime(){}
void start(int world_rank);
void  start();
void end(int world_rank);
void end();
};
void mytime::start(int world_rank)
{	t1=time(NULL);}

void mytime::start()
{	t1=time(NULL);

}

void mytime::end(int world_rank)
{  if(world_rank==0){
    t2=time(NULL);                 
    printf ("time = %d secs \n =%d hours+%d mins+%d secs\n", (t2 - t1),int((t2-t1)/3600),int((t2-t1)%3600/60),int((t2-t1)%3600%60)); 
}
}

void mytime::end()
{
    t2=time(NULL);                 
    printf ("time = %d secs \n =%d hours+%d mins+%d secs\n", (t2 - t1),int((t2-t1)/3600),int((t2-t1)%3600/60),int((t2-t1)%3600%60)); 
}
