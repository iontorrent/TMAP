// Coder: ACRush
// Submission: 27
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=19849563&subnum=27
//#define LOCAL
//#define DEBUG
//#define PRINT

#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <string>
#include <cstring>
#include <ctime>

//#ifndef LOCAL
#include <emmintrin.h>
#include <pmmintrin.h>
#include <xmmintrin.h>
//#endif
#include "solution5.h"

using namespace std;

#define SIZE(A) ((int)A.size())
#define LENGTH(A) ((int)A.length())
typedef long long int64;
typedef unsigned long long uint64;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef pair<int,int> ipair;
#define MP(A,B) make_pair(A,B)
typedef vector<string> VS;
typedef vector<int> VI;
#define two(X) (1<<(X))
#define contain(S,X) ((S&two(X))>0)
const double pi=acos(-1.0);
inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }
static int opt,query_end,target_end,n_best;

int LOG(const char* s,...)
{
#ifdef PRINT
  va_list va;
  va_start(va,s);
  FILE *f=fopen("test.log","a");
  int ret=vfprintf(f,s,va);
  fclose(f);
  va_end(va);
  return ret;
#else
  return 0;
#endif
}

#ifdef DEBUG
#define ASSERT(_Expression) (void)((!!(_Expression)||(__assert((#_Expression),__LINE__),0)));
void __assert(const char *_Message,const unsigned _Line)
{
  printf("ASSERTION FAILED\n");
  printf("Message = %s\n",_Message);
  printf("Line = %u\n",_Line);
  exit(0);
}
#endif

#ifdef LOCAL
#ifdef _MSC_VER
#define TIMES_PER_SEC (1.75e9)
#else
#define TIMES_PER_SEC (1.6508566e9)
#endif
#else
#define TIMES_PER_SEC (3.589e9)
#endif

#ifdef LOCAL
class Timer
{
public:
  static uint64 timeUsed[1000];
  int id;
  uint64 startTime;
  Timer(int id=0)
    {
      this->id=id;
      startTime=rdtsc();
    }
  ~Timer()
    {
      timeUsed[id]+=(rdtsc()-startTime);
    }
  static void show()
    {
      for (int i=0;i<1000;i++) if (timeUsed[i]>0)
        {
          char str[100];
          sprintf(str,"%.6lf",timeUsed[i]/TIMES_PER_SEC);
          string s=str;
          if (LENGTH(s)<15) s=" "+s;
          printf("%4d %s\n",i,s.c_str());
        }
    }
};
uint64 Timer::timeUsed[1000]={0};

class Counter
{
public:
  static uint64 cnt[1000];
  Counter(int id=0)
    {
      cnt[id]++;
    }
  ~Counter()
    {
    }
  static void show()
    {
      for (int i=0;i<1000;i++) if (cnt[i]>0)
        printf("Counter %d %d\n",i,cnt[i]);
    }
};
uint64 Counter::cnt[1000]={0};
#endif

const int MAX_LENGTH=1040;
const int INF=1073741824;

inline int max(int a,int b)
{
  return (a>b)?a:b;
}
inline int min(int a,int b)
{
  return (a<b)?a:b;
}
inline int max(int a,int b,int c)
{
  return (a>b)?((a>c)?a:c):((b>c)?b:c);
}

#ifdef LOCAL
int __opt;
#endif

#ifdef LOCAL
#define ALIGN16 __declspec(align(16))
#endif

static int m,n;
char target[MAX_LENGTH],query[MAX_LENGTH];
int a[MAX_LENGTH],b[MAX_LENGTH],ca[MAX_LENGTH],cb[MAX_LENGTH];
int task_id;
int match_score,mismatch_score,gap_open,gap_extension;
int direction;
uint batch8[1024];
uint mismatch_score_batch16;

struct dpstate
{
  int idx,f,h;
};
int sizeq[2];
dpstate dp_buffer[MAX_LENGTH*3+100];

int q[3][MAX_LENGTH];
struct dpstate_large
{
  int f,h;
};
dpstate_large dp_large[MAX_LENGTH][MAX_LENGTH];
int counter_task2;
bool vf[MAX_LENGTH][MAX_LENGTH];
bool vh[MAX_LENGTH][MAX_LENGTH];

#ifdef LOCAL
ALIGN16 short hbuffer16[MAX_LENGTH*3+100];
ALIGN16 uchar hbuffer8[MAX_LENGTH*3+100];
#else
short hbuffer16[MAX_LENGTH*3+100] __attribute__ ((aligned (16)));
uchar hbuffer8[MAX_LENGTH*3+100] __attribute__ ((aligned (16)));
#endif

#ifdef LOCAL
ALIGN16 int int_buffer[4];
ALIGN16 short short_buffer[8];
ALIGN16 uchar uchar_buffer[16];
#else
int int_buffer[4] __attribute__ ((aligned (16)));
short short_buffer[8] __attribute__ ((aligned (16)));
uchar uchar_buffer[16] __attribute__ ((aligned (16)));
#endif


int vdeg[256],vgraph[256][MAX_LENGTH];
int vfirst[256],vnext[MAX_LENGTH];

int testcase_id=0;
int measure_id;

int est_opt;
int first_j,last_j;
char result[1024];

bool valid0[256];
int encode_idx[128];

void init()
{
  encode_idx['A']=0;
  encode_idx['C']=1;
  encode_idx['G']=2;
  encode_idx['T']=3;
  counter_task2=0;
  for (int set=0;set<256;set++)
    {
      valid0[set]=true;
      int cnt=0;
      if ((set&3)==0) cnt++;
      if ((set&(3<<2))==0) cnt++;
      if ((set&(3<<4))==0) cnt++;
      if ((set&(3<<6))==0) cnt++;
      if (cnt<=2) valid0[set]=false;
    }
}

void solve();

int get_task_id(int qsc,int qec)
{
  if (qsc==1 && qec==1) return 1;
  if (qsc==0 && qec==1) return 2;
  if (qsc==1 && qec==0) return 3;
  if (qsc==0 && qec==0) return 4;
  return 0;
}

Solution5::Solution5() {
}

#ifdef LOCAL
int process(string& s_target, string& s_query, int qsc, int qec, int v_match_score, 
               int v_mismatch_score, int v_gap_open, int v_gap_extension, int v_direction,
               int *_opt, int *_te, int *_qe, int *_n_best)
#else
int Solution5::
process(string& s_target, string& s_query, int qsc, int qec, int v_match_score, 
        int v_mismatch_score, int v_gap_open, int v_gap_extension, int v_direction,
        int *_opt, int *_te, int *_qe, int *_n_best)
#endif
{
#ifdef LOCAL
  Timer t99(99);
#endif
  if (testcase_id==0)
    {
      init();
      match_score=v_match_score;
      mismatch_score=v_mismatch_score;
      gap_open=v_gap_open;
      gap_extension=v_gap_extension;
      if (match_score==1 && mismatch_score==-3 && gap_open==-5 && gap_extension==-2) measure_id=0;
      if (match_score==4 && mismatch_score==-5 && gap_open==-8 && gap_extension==-2) measure_id=1;
      uchar u[4];
      for (int set=0;set<1024;set++)
        {
          int last=set&3;
          for (int i=0;i<4;i++) u[i]=(((set>>(i+i+2))&3)==last)?match_score:mismatch_score;
          memcpy(&batch8[set],u,sizeof(uint));
        }
      short s[2];
      for (int i=0;i<2;i++) s[i]=mismatch_score;
      memcpy(&mismatch_score_batch16,s,sizeof(uint));
    }
  int v_task_idx=get_task_id(qsc,qec);
  int l_query=LENGTH(s_query);
  int l_target=LENGTH(s_target);
  const char *v_query=s_query.c_str();
  const char *v_target=s_target.c_str();
  est_opt=0;
  if (testcase_id>0 && 
      direction==v_direction && 
      task_id==v_task_idx &&
      n==l_query && 
      m>=l_target &&
      memcmp(v_query,query,n*sizeof(char))==0)
    {
#ifdef LOCAL
      Timer t0(0);
#endif
      int mid_l=m-l_target;
      int max_l=mid_l+l_target/10;
      int ck_l=l_target;
      for (int src=0;src<=max_l;src++)
        {
          if (src>mid_l) ck_l--;
          if (memcmp(v_target,target+src,ck_l*sizeof(char))==0)
            {
              int dest=src+l_target-1;
              //sscanf(result,"%d%d%d%d",&opt,&query_end,&target_end,&n_best);
              if (dest<m && first_j>=src && last_j<=dest)
                {
#ifdef LOCAL
                  Counter c91(91);
#endif
                  char new_result[1024];
                  //sprintf(new_result,"%d %d %d %d",opt,query_end,target_end-src,n_best);
                  target_end -= src;
                  testcase_id++;
                  (*_opt) = opt;
                  (*_qe) = query_end;
                  (*_te) = target_end;
                  (*_n_best) = n_best;
                  return opt;
                }
              if (dest>m-1) dest=m-1;
              if (first_j<=dest && src<=last_j)
                {
#ifdef LOCAL
                  Counter c92(92);
#endif
                  int range=last_j-first_j+1;
                  int overlap=min(last_j,dest)-max(first_j,src)+1;
                  if (task_id==1)
                    {
                      est_opt=max(0,opt-(range-overlap)*match_score);
                    }
                  else if (task_id==2 || task_id==3)
                    {
                      if (overlap<=0) 
                        est_opt=0;
                      else
                        {
                          int t=opt-(range-overlap)*match_score;
                          if (task_id==2 && first_j<src) t+=gap_open+gap_extension*(src-first_j);
                          if (task_id==3 && last_j>dest) t+=gap_open+gap_extension*(last_j-dest);
                          est_opt=max(0,t);
                        }
                    }
                }
              break;
            }
        }
    }
#ifdef LOCAL
  //est_opt=__opt;
#endif
  n=l_query;
  memcpy(query,v_query,n*sizeof(char));
  m=l_target;
  memcpy(target,v_target,m*sizeof(char));
  task_id=v_task_idx;
  direction=v_direction;
  solve();
  testcase_id++;
  //return result;
  (*_opt) = opt;
  (*_qe) = query_end;
  (*_te) = target_end;
  (*_n_best) = n_best;
  return opt;
}

/*
   inline int encode(char c)
   {
   switch (c)
   {
   case 'A': return 0;
   case 'C': return 1;
   case 'G': return 2;
   case 'T': return 3;
   }
   return -1;
   }
   */
#define encode(c) encode_idx[c]

void shift_right16(__m128i &a,short first)
{
  __m128i *p=(__m128i*)short_buffer;
  _mm_store_si128(p,a);
  for (int i=7;i>0;i--) short_buffer[i]=short_buffer[i-1];
  short_buffer[0]=first;
  a=_mm_load_si128(p);
}

void shift_right8(__m128i &a,uchar first)
{
  __m128i *p=(__m128i*)uchar_buffer;
  _mm_store_si128(p,a);
  for (int i=15;i>0;i--) uchar_buffer[i]=uchar_buffer[i-1];
  uchar_buffer[0]=first;
  a=_mm_load_si128(p);
}

int shuffle[MAX_LENGTH];
#ifdef LOCAL
ALIGN16 uchar cost8[4][MAX_LENGTH];
ALIGN16 short cost16[4][MAX_LENGTH];
#else
uchar cost8[4][MAX_LENGTH] __attribute__ ((aligned (16)));
short cost16[4][MAX_LENGTH] __attribute__ ((aligned (16)));
#endif

short transfer_h[MAX_LENGTH],transfer_e[MAX_LENGTH];

#ifdef LOCAL
int cc1=0,cc2=0,cc3=0;
#endif

#define ZERO (45)
#define LBOUND1 (30)
#define LBOUND2 (15)

void brute_force()
{
#ifdef LOCAL
  Timer t95(95);
#endif
    {
#ifdef LOCAL
      Timer t70(70);
#endif
      for (int i=0;i<n;i++) a[i]=encode(query[i]);
      for (int j=0;j<m;j++) b[j]=encode(target[j]);
    }
  gap_open+=gap_extension;
  //int opt=1,min_last_j,max_last_j,n_best,query_end,target_end;
  int min_last_j,max_last_j; 
  opt=-1;
  if (est_opt-1>opt) opt=est_opt-1;
  int start_i=0;
  int m_ext_g=60;
  if (task_id==1 || task_id==3)
    {
      int mk=(m+15)/16;
      int rounds=m_ext_g/mk+1;
      int length=(mk<<4);
      for (int pos=0,src=0;src<mk;src++) for (int j=src;j<length;j+=mk) shuffle[pos++]=j;
      uint *p0=(uint*)cost8[0];
      uint *p1=(uint*)cost8[1];
      uint *p2=(uint*)cost8[2];
      uint *p3=(uint*)cost8[3];
      for (int j=m;j<length;j++) b[j]=0;
      for (int j=0;j<length;j+=4)
        {
          int set=(b[shuffle[j]]<<2)|(b[shuffle[j+1]]<<4)|(b[shuffle[j+2]]<<6)|(b[shuffle[j+3]]<<8);
          *p0=batch8[set]; *p1=batch8[set|1]; *p2=batch8[set|2]; *p3=batch8[set|3];
          p0++; p1++; p2++; p3++;
        }
      char c_mismatch_score=mismatch_score;
      for (int pos=15;pos<length;pos+=16) if (shuffle[pos]>=m) cost8[0][pos]=c_mismatch_score;
      __m128i m_gap_open=_mm_set1_epi8(gap_open);
      __m128i m_gap_extension=_mm_set1_epi8(gap_extension);
      __m128i m_gap_extension_x2=_mm_set1_epi8(gap_extension*2);
      __m128i m_gap_extension_x3=_mm_set1_epi8(gap_extension*3);
      __m128i m_zero=_mm_set1_epi8(ZERO);
      __m128i m_lbound1=_mm_set1_epi8(LBOUND1);
      __m128i m_lbound2=_mm_set1_epi8(LBOUND2);
      __m128i *p_buffer=(__m128i*)uchar_buffer;
      uchar *h[2];
      h[0]=hbuffer8;
      h[1]=hbuffer8+length;
      uchar *e=hbuffer8+length+length;
      for (int j=0;j<length;j++) { h[0][j]=ZERO; e[j]=ZERO+gap_open; }
      __m128i *vE=(__m128i*)e;
      __m128i *ck_zero=(__m128i*)int_buffer;
      __m128i *vHLoad=(__m128i*)h[1];
      __m128i *vHStore=(__m128i*)h[0];
      bool overflow=false;
      for (int i=0;i<n;i++)
        {
          __m128i *vProfile=(__m128i*)cost8[a[i]];
          __m128i vF=m_lbound1;
          __m128i m01,m02,m03,m04,m05,m06,m07,m08,m09;
          __m128i m11,m12,m13,m14,m15,m16,m17,m18,m19;
          __m128i vH2,vF2;
          __m128i *tmp=vHLoad;
          vHLoad=vHStore;
          vHStore=tmp;
          m01=_mm_load_si128(vHLoad+mk-1);
          __m128i vH=_mm_max_epu8(m01,m_zero);
          shift_right8(vH,ZERO);
          int j=0;
          for (;j+1<mk;j+=2)            
            {
              m01=_mm_load_si128(vProfile+j);
              m02=_mm_load_si128(vE+j);
              m09=_mm_load_si128(vHLoad+j);
              m03=_mm_add_epi8(vH,m01);
              m04=_mm_max_epu8(m02,m03);
              vH2=_mm_max_epu8(m09,m_zero);
              vH=_mm_max_epu8(vF,m04);
              m08=_mm_add_epi8(vF,m_gap_extension);
              m05=_mm_add_epi8(vH,m_gap_open);
              vF2=_mm_max_epu8(m05,m08);

              m11=_mm_load_si128(vProfile+j+1);
              m12=_mm_load_si128(vE+j+1);
              m13=_mm_add_epi8(vH2,m11);
              m14=_mm_max_epu8(m12,m13);
              vH2=_mm_max_epu8(vF2,m14);

              _mm_store_si128(vHStore+j,vH);
              m19=_mm_load_si128(vHLoad+j+1);
              vH=_mm_max_epu8(m19,m_zero);
              m15=_mm_add_epi8(vH2,m_gap_open);
              m18=_mm_add_epi8(vF2,m_gap_extension);
              m16=_mm_add_epi8(m12,m_gap_extension);
              m06=_mm_add_epi8(m02,m_gap_extension);
              vF=_mm_max_epu8(m15,m18);
              m17=_mm_max_epu8(m15,m16);
              m07=_mm_max_epu8(m05,m06);

              _mm_store_si128(vE+j,m07);
              _mm_store_si128(vHStore+j+1,vH2);
              _mm_store_si128(vE+j+1,m17);
            }
          if (j<mk)
            {
              m01=_mm_load_si128(vProfile+j);
              m02=_mm_load_si128(vE+j);
              m06=_mm_add_epi8(m02,m_gap_extension);
              m08=_mm_add_epi8(vF,m_gap_extension);
              m03=_mm_add_epi8(vH,m01);
              m04=_mm_max_epu8(m02,m03);
              vH=_mm_max_epu8(vF,m04);
              m05=_mm_add_epi8(vH,m_gap_open);
              vF=_mm_max_epu8(m05,m08);
              m07=_mm_max_epu8(m05,m06);
              _mm_store_si128(vHStore+j,vH);
              _mm_store_si128(vE+j,m07);
              m09=_mm_load_si128(vHLoad+j);
              vH=_mm_max_epu8(m09,m_zero);
            }
          shift_right8(vF,LBOUND2);
          for (int counter=0;counter<rounds;counter++)
            {
              bool flag=false;
              for (j=0;j+3<mk;j+=4)
                if (j&15)
                  {
                    m01=_mm_load_si128(vHStore+j);
                    m12=_mm_add_epi8(vF,m_gap_extension);
                    m06=_mm_add_epi8(vF,m_gap_extension_x2);
                    m16=_mm_add_epi8(vF,m_gap_extension_x3);
                    m11=_mm_load_si128(vHStore+j+1);
                    m07=_mm_load_si128(vHStore+j+2);
                    m17=_mm_load_si128(vHStore+j+3);
                    m04=_mm_max_epu8(m01,vF);
                    m14=_mm_max_epu8(m11,m12);
                    m08=_mm_max_epu8(m06,m07);
                    m18=_mm_max_epu8(m16,m17);
                    _mm_store_si128(vHStore+j,m04);
                    _mm_store_si128(vHStore+j+1,m14);
                    _mm_store_si128(vHStore+j+2,m08);
                    _mm_store_si128(vHStore+j+3,m18);
                    m09=_mm_add_epi8(m16,m_gap_extension);
                    vF=_mm_max_epu8(m09,m_lbound2);
                  }
                else
                  {
                    m01=_mm_load_si128(vHStore+j);
                    m02=_mm_add_epi8(m01,m_gap_open);
                    m03=_mm_max_epu8(vF,m02);
                    m06=_mm_cmpeq_epi8(vF,m03);
                    _mm_store_si128(ck_zero,m06);
                    if (int_buffer[0]==0 && int_buffer[1]==0 && int_buffer[2]==0 && int_buffer[3]==0) { flag=true; break; }
                    m12=_mm_add_epi8(vF,m_gap_extension);
                    m06=_mm_add_epi8(vF,m_gap_extension_x2);
                    m16=_mm_add_epi8(vF,m_gap_extension_x3);
                    m11=_mm_load_si128(vHStore+j+1);
                    m07=_mm_load_si128(vHStore+j+2);
                    m17=_mm_load_si128(vHStore+j+3);
                    m04=_mm_max_epu8(m01,vF);
                    m14=_mm_max_epu8(m11,m12);
                    m08=_mm_max_epu8(m06,m07);
                    m18=_mm_max_epu8(m16,m17);
                    _mm_store_si128(vHStore+j,m04);
                    _mm_store_si128(vHStore+j+1,m14);
                    _mm_store_si128(vHStore+j+2,m08);
                    _mm_store_si128(vHStore+j+3,m18);
                    m09=_mm_add_epi8(m16,m_gap_extension);
                    vF=_mm_max_epu8(m09,m_lbound2);
                  }
              if (flag) break;
              for (;j<mk;j++)
                {
                  m01=_mm_load_si128(vHStore+j);
                  if (j+1==mk)
                    {
                      m02=_mm_add_epi8(m01,m_gap_open);
                      m03=_mm_max_epu8(vF,m02);
                      m06=_mm_cmpeq_epi8(vF,m03);
                      _mm_store_si128(ck_zero,m06);
                      if (int_buffer[0]==0 && int_buffer[1]==0 && int_buffer[2]==0 && int_buffer[3]==0) { flag=true; break; }
                    }
                  m04=_mm_max_epu8(m01,vF);
                  _mm_store_si128(vHStore+j,m04);
                  m05=_mm_add_epi8(vF,m_gap_extension);
                  vF=_mm_max_epu8(m05,m_lbound2);
                }
              if (flag) break;
              shift_right8(vF,LBOUND2);
            }
          m01=m_zero; m02=m_zero; m03=m_zero; m04=m_zero; 
          m05=m_zero; m06=m_zero; m07=m_zero; m08=m_zero; 
          for (j=0;j+7<mk;j+=8)
            {
              m11=_mm_load_si128(vHStore+j); m12=_mm_load_si128(vHStore+j+1);
              m13=_mm_load_si128(vHStore+j+2); m14=_mm_load_si128(vHStore+j+3);
              m15=_mm_load_si128(vHStore+j+4); m16=_mm_load_si128(vHStore+j+5);
              m17=_mm_load_si128(vHStore+j+6); m18=_mm_load_si128(vHStore+j+7);
              m01=_mm_max_epu8(m01,m11); m02=_mm_max_epu8(m02,m12);
              m03=_mm_max_epu8(m03,m13); m04=_mm_max_epu8(m04,m14);
              m05=_mm_max_epu8(m05,m15); m06=_mm_max_epu8(m06,m16);
              m07=_mm_max_epu8(m07,m17); m08=_mm_max_epu8(m08,m18);
            }
          for (;j<mk;j++)
            {
              m11=_mm_load_si128(vHStore+j);
              m01=_mm_max_epu8(m01,m11);
            }
          m11=_mm_max_epu8(m01,m02);
          m13=_mm_max_epu8(m03,m04);
          m15=_mm_max_epu8(m05,m06);
          m17=_mm_max_epu8(m07,m08);
          m12=_mm_max_epu8(m11,m13);
          m14=_mm_max_epu8(m15,m17);
          m16=_mm_max_epu8(m12,m14);
          _mm_store_si128(p_buffer,m16);
          int local_opt=0;
          for (int d=0;d<16;d++) if ((int)uchar_buffer[d]-ZERO>local_opt) local_opt=uchar_buffer[d]-ZERO;
          if (task_id==1 && local_opt>=opt)
            {
              uchar *h0=h[(i+1)&1];
              for (int d=0;d<16;d++) if ((int)uchar_buffer[d]-ZERO==local_opt)
                for (int pos_j=shuffle[d],j=d;pos_j<m && j<length;j+=16,pos_j++)
                  {
                    int t=(int)h0[j]-ZERO;
                    if (t>=opt) 
                      if (t>opt) 
                        opt=t,n_best=1,query_end=i,target_end=min_last_j=max_last_j=pos_j;
                      else 
                        {
                          n_best++;
                          if (direction==1) query_end=i,target_end=pos_j;
                          if (pos_j>max_last_j) max_last_j=pos_j;
                          if (pos_j<min_last_j) min_last_j=pos_j;
                        }
                  }
            }
          if (local_opt>190) 
            { 
              overflow=true; 
              start_i=i+1; 
              uchar *h0=h[(i+1)&1];
              for (int j=0;j<length;j++) 
                {
                  int pos=shuffle[j];
                  transfer_h[shuffle[j]]=h0[j]-ZERO;
                  transfer_e[shuffle[j]]=e[j]-ZERO;
                }
              break; 
            }
        }
      if (!overflow)
        {
          if (task_id==3)
            {
              query_end=n-1;
              for (int j=0;j<length;j++) if (shuffle[j]<m)
                {
                  int t=h[n&1][j]-ZERO;
                  if (t>=opt) 
                    if (t>opt) 
                      opt=t,n_best=1,min_last_j=max_last_j=shuffle[j];
                    else 
                      {
                        n_best++;
                        if (shuffle[j]>max_last_j) max_last_j=shuffle[j];
                        if (shuffle[j]<min_last_j) min_last_j=shuffle[j];
                      }
                }
              target_end=(direction==1)?max_last_j:min_last_j;
            }
          first_j=min_last_j-n-(n*match_score-opt)/(-gap_extension)-2;
          if (first_j<0) first_j=0;
          last_j=max_last_j;
          //sprintf(result,"%d %d %d %d",opt,query_end,target_end,n_best);
          gap_open-=gap_extension;
          return;
        }
    }
  int mk=(m+7)/8;
  int rounds=m_ext_g/mk+1;
  int length=(mk<<3);
    {
#ifdef LOCAL
      Timer t71(71);
#endif
      for (int c=0;c<4;c++)
        {
          uint *tmp=(uint*)cost16[c];
          for (int j=(length>>1)-1;j>=0;j--) tmp[j]=mismatch_score_batch16;
        }
      short s_match_score=match_score;
      for (int pos=0,src=0;src<mk;src++) for (int j=src;j<length;j+=mk) 
        {
          shuffle[pos]=j;
          if (j<m) cost16[b[j]][pos]=s_match_score;
          pos++;
        }
    }
  __m128i m_gap_open=_mm_set1_epi16(gap_open);
  __m128i m_gap_extension=_mm_set1_epi16(gap_extension);
  __m128i m_gap_extension_x2=_mm_set1_epi16(gap_extension*2);
  __m128i m_gap_extension_x3=_mm_set1_epi16(gap_extension*3);
  __m128i m_zero=_mm_set1_epi16(0);
  short *h[2];
  h[0]=hbuffer16;
  h[1]=hbuffer16+length;
  short *e=hbuffer16+length+length;
    {
#ifdef LOCAL
      Timer t72(72);
#endif
      if (start_i==0)
        for (int j=0;j<length;j++) { h[0][j]=0; e[j]=gap_open; }
      else
        {
          short *nexth=h[start_i&1];
          for (int j=0;j<length;j++) { nexth[j]=transfer_h[shuffle[j]]; e[j]=transfer_e[shuffle[j]]; }
        }
    }
  __m128i *vE=(__m128i*)e;
  __m128i *ck_zero=(__m128i*)int_buffer;
  __m128i *p_buffer=(__m128i*)short_buffer;
  __m128i *vHLoad=(__m128i*)h[(start_i+1)&1];
  __m128i *vHStore=(__m128i*)h[start_i&1];
    {
#ifdef LOCAL
      Timer t73(73);
#endif
      for (int i=start_i;i<n;i++)
        {
          __m128i *vProfile=(__m128i*)cost16[a[i]];
          __m128i vF=_mm_set1_epi16(-10000);
          __m128i m01,m02,m03,m04,m05,m06,m07,m08,m09;
          __m128i m11,m12,m13,m14,m15,m16,m17,m18,m19;
          __m128i vH2,vF2;
          __m128i *tmp=vHLoad;
          vHLoad=vHStore;
          vHStore=tmp;
          if (task_id==1 || task_id==3)
            {
              m01=_mm_load_si128(vHLoad+mk-1);
              __m128i vH=_mm_max_epi16(m01,m_zero);
              shift_right16(vH,0);
              int j=0;
              for (;j+1<mk;j+=2)
                {
                  m01=_mm_load_si128(vProfile+j);
                  m02=_mm_load_si128(vE+j);
                  m09=_mm_load_si128(vHLoad+j);
                  m03=_mm_add_epi16(vH,m01);
                  m04=_mm_max_epi16(m02,m03);
                  vH2=_mm_max_epi16(m09,m_zero);
                  vH=_mm_max_epi16(vF,m04);
                  m08=_mm_add_epi16(vF,m_gap_extension);
                  m05=_mm_add_epi16(vH,m_gap_open);
                  vF2=_mm_max_epi16(m05,m08);

                  m11=_mm_load_si128(vProfile+j+1);
                  m12=_mm_load_si128(vE+j+1);
                  m13=_mm_add_epi16(vH2,m11);
                  m14=_mm_max_epi16(m12,m13);
                  vH2=_mm_max_epi16(vF2,m14);

                  _mm_store_si128(vHStore+j,vH);
                  m19=_mm_load_si128(vHLoad+j+1);
                  vH=_mm_max_epi16(m19,m_zero);
                  m15=_mm_add_epi16(vH2,m_gap_open);
                  m18=_mm_add_epi16(vF2,m_gap_extension);
                  m16=_mm_add_epi16(m12,m_gap_extension);
                  m06=_mm_add_epi16(m02,m_gap_extension);
                  vF=_mm_max_epi16(m15,m18);
                  m17=_mm_max_epi16(m15,m16);
                  m07=_mm_max_epi16(m05,m06);

                  _mm_store_si128(vE+j,m07);
                  _mm_store_si128(vHStore+j+1,vH2);
                  _mm_store_si128(vE+j+1,m17);
                }
              if (j<mk)
                {
                  m01=_mm_load_si128(vProfile+j);
                  m02=_mm_load_si128(vE+j);
                  m06=_mm_add_epi16(m02,m_gap_extension);
                  m08=_mm_add_epi16(vF,m_gap_extension);
                  m03=_mm_add_epi16(vH,m01);
                  m04=_mm_max_epi16(m02,m03);
                  vH=_mm_max_epi16(vF,m04);
                  m05=_mm_add_epi16(vH,m_gap_open);
                  vF=_mm_max_epi16(m05,m08);
                  m07=_mm_max_epi16(m05,m06);
                  _mm_store_si128(vHStore+j,vH);
                  _mm_store_si128(vE+j,m07);
                  m09=_mm_load_si128(vHLoad+j);
                  vH=_mm_max_epi16(m09,m_zero);
                }
            }
          else
            {
              __m128i vH=_mm_load_si128(vHLoad+mk-1);
              shift_right16(vH,(i==0)?0:(gap_open+(i-1)*gap_extension));
              int j=0;
              for (;j+1<mk;j+=2)
                {
                  m01=_mm_load_si128(vProfile+j);
                  m02=_mm_load_si128(vE+j);
                  m03=_mm_add_epi16(vH,m01);
                  m04=_mm_max_epi16(m02,m03);
                  vH=_mm_max_epi16(vF,m04);
                  m08=_mm_add_epi16(vF,m_gap_extension);
                  m05=_mm_add_epi16(vH,m_gap_open);
                  vF2=_mm_max_epi16(m05,m08);

                  vH2=_mm_load_si128(vHLoad+j);
                  m11=_mm_load_si128(vProfile+j+1);
                  m12=_mm_load_si128(vE+j+1);
                  m13=_mm_add_epi16(vH2,m11);
                  m14=_mm_max_epi16(m12,m13);
                  vH2=_mm_max_epi16(vF2,m14);

                  _mm_store_si128(vHStore+j,vH);
                  m15=_mm_add_epi16(vH2,m_gap_open);
                  m18=_mm_add_epi16(vF2,m_gap_extension);
                  m16=_mm_add_epi16(m12,m_gap_extension);
                  m06=_mm_add_epi16(m02,m_gap_extension);
                  vF=_mm_max_epi16(m15,m18);
                  m17=_mm_max_epi16(m15,m16);
                  m07=_mm_max_epi16(m05,m06);
                  vH=_mm_load_si128(vHLoad+j+1);

                  _mm_store_si128(vE+j,m07);
                  _mm_store_si128(vHStore+j+1,vH2);
                  _mm_store_si128(vE+j+1,m17);
                }
              if (j<mk)
                {
                  m01=_mm_load_si128(vProfile+j);
                  m02=_mm_load_si128(vE+j);
                  m06=_mm_add_epi16(m02,m_gap_extension);
                  m08=_mm_add_epi16(vF,m_gap_extension);
                  m03=_mm_add_epi16(vH,m01);
                  m04=_mm_max_epi16(m02,m03);
                  vH=_mm_max_epi16(vF,m04);
                  m05=_mm_add_epi16(vH,m_gap_open);
                  vF=_mm_max_epi16(m05,m08);
                  m07=_mm_max_epi16(m05,m06);
                  _mm_store_si128(vHStore+j,vH);
                  _mm_store_si128(vE+j,m07);
                  vH=_mm_load_si128(vHLoad+j);
                }
            }
          if (task_id==4)
            {
              shift_right16(vF,-10000);
              for (int counter=0;counter<rounds;counter++)
                {
                  int j=0;
                  for (;j+3<mk;j+=4)
                    {
                      m01=_mm_load_si128(vHStore+j);
                      m13=_mm_add_epi16(vF,m_gap_extension);
                      m04=_mm_add_epi16(vF,m_gap_extension_x2);
                      m14=_mm_add_epi16(vF,m_gap_extension_x3);
                      m11=_mm_load_si128(vHStore+j+1);
                      m05=_mm_load_si128(vHStore+j+2);
                      m15=_mm_load_si128(vHStore+j+3);
                      m02=_mm_max_epi16(m01,vF);
                      m12=_mm_max_epi16(m11,m13);
                      m06=_mm_max_epi16(m05,m04);
                      m16=_mm_max_epi16(m15,m14);
                      _mm_store_si128(vHStore+j,m02);
                      _mm_store_si128(vHStore+j+1,m12);
                      _mm_store_si128(vHStore+j+2,m06);
                      _mm_store_si128(vHStore+j+3,m16);
                      vF=_mm_add_epi16(m14,m_gap_extension);
                    }
                  for (;j<mk;j++)
                    {
                      m01=_mm_load_si128(vHStore+j);
                      m04=_mm_max_epi16(m01,vF);
                      _mm_store_si128(vHStore+j,m04);
                      vF=_mm_add_epi16(vF,m_gap_extension);
                    }
                  shift_right16(vF,-10000);
                }
            }
          else
            {
              shift_right16(vF,-10000);
              for (int counter=0;counter<rounds;counter++)
                {
                  bool flag=false;
                  int j=0;
                  for (;j+3<mk;j+=4)
                    if (j&15)
                      {
                        m01=_mm_load_si128(vHStore+j);
                        m12=_mm_add_epi16(vF,m_gap_extension);
                        m06=_mm_add_epi16(vF,m_gap_extension_x2);
                        m16=_mm_add_epi16(vF,m_gap_extension_x3);
                        m11=_mm_load_si128(vHStore+j+1);
                        m07=_mm_load_si128(vHStore+j+2);
                        m17=_mm_load_si128(vHStore+j+3);
                        m04=_mm_max_epi16(m01,vF);
                        m14=_mm_max_epi16(m11,m12);
                        m08=_mm_max_epi16(m06,m07);
                        m18=_mm_max_epi16(m16,m17);
                        _mm_store_si128(vHStore+j,m04);
                        _mm_store_si128(vHStore+j+1,m14);
                        _mm_store_si128(vHStore+j+2,m08);
                        _mm_store_si128(vHStore+j+3,m18);
                        vF=_mm_add_epi16(m16,m_gap_extension);
                      }
                    else
                      {
                        m01=_mm_load_si128(vHStore+j);
                        m02=_mm_add_epi16(m01,m_gap_open);
                        m03=_mm_cmpgt_epi16(vF,m02);
                        _mm_store_si128(ck_zero,m03);
                        if (int_buffer[0]==0 && int_buffer[1]==0 && int_buffer[2]==0 && int_buffer[3]==0) { flag=true; break; }
                        m12=_mm_add_epi16(vF,m_gap_extension);
                        m06=_mm_add_epi16(vF,m_gap_extension_x2);
                        m16=_mm_add_epi16(vF,m_gap_extension_x3);
                        m11=_mm_load_si128(vHStore+j+1);
                        m07=_mm_load_si128(vHStore+j+2);
                        m17=_mm_load_si128(vHStore+j+3);
                        m04=_mm_max_epi16(m01,vF);
                        m14=_mm_max_epi16(m11,m12);
                        m08=_mm_max_epi16(m06,m07);
                        m18=_mm_max_epi16(m16,m17);
                        _mm_store_si128(vHStore+j,m04);
                        _mm_store_si128(vHStore+j+1,m14);
                        _mm_store_si128(vHStore+j+2,m08);
                        _mm_store_si128(vHStore+j+3,m18);
                        vF=_mm_add_epi16(m16,m_gap_extension);
                      }
                  if (flag) break;
                  for (;j<mk;j++)
                    {
                      m01=_mm_load_si128(vHStore+j);
                      if (j+1==mk)
                        {
                          m02=_mm_add_epi16(m01,m_gap_open);
                          m03=_mm_cmpgt_epi16(vF,m02);
                          _mm_store_si128(ck_zero,m03);
                          if (int_buffer[0]==0 && int_buffer[1]==0 && int_buffer[2]==0 && int_buffer[3]==0) { flag=true; break; }
                        }
                      m04=_mm_max_epi16(m01,vF);
                      _mm_store_si128(vHStore+j,m04);
                      vF=_mm_add_epi16(vF,m_gap_extension);
                    }
                  if (flag) break;
                  shift_right16(vF,-10000);
                }
            }
          if (task_id==1 || task_id==2)
            {
#ifdef LOCAL
              Timer t74(74);
#endif
              m01=m_zero; m02=m_zero; m03=m_zero; m04=m_zero; 
              m05=m_zero; m06=m_zero; m07=m_zero; m08=m_zero; 
              int j=0;    
              for (;j+7<mk;j+=8)
                {
                  m11=_mm_load_si128(vHStore+j); m12=_mm_load_si128(vHStore+j+1);
                  m13=_mm_load_si128(vHStore+j+2); m14=_mm_load_si128(vHStore+j+3);
                  m15=_mm_load_si128(vHStore+j+4); m16=_mm_load_si128(vHStore+j+5);
                  m17=_mm_load_si128(vHStore+j+6); m18=_mm_load_si128(vHStore+j+7);
                  m01=_mm_max_epi16(m01,m11); m02=_mm_max_epi16(m02,m12);
                  m03=_mm_max_epi16(m03,m13); m04=_mm_max_epi16(m04,m14);
                  m05=_mm_max_epi16(m05,m15); m06=_mm_max_epi16(m06,m16);
                  m07=_mm_max_epi16(m07,m17); m08=_mm_max_epi16(m08,m18);
                }
              for (;j<mk;j++)
                {
                  m11=_mm_load_si128(vHStore+j);
                  m01=_mm_max_epi16(m01,m11);
                }
              m11=_mm_max_epi16(m01,m02);
              m13=_mm_max_epi16(m03,m04);
              m15=_mm_max_epi16(m05,m06);
              m17=_mm_max_epi16(m07,m08);
              m12=_mm_max_epi16(m11,m13);
              m14=_mm_max_epi16(m15,m17);
              m16=_mm_max_epi16(m12,m14);
              _mm_store_si128(p_buffer,m16);
              int local_opt=0;
              for (int d=0;d<8;d++) if ((int)short_buffer[d]>local_opt) local_opt=short_buffer[d];
              if (local_opt>0 && local_opt>=opt)
                {
                  short *h0=h[(i+1)&1];
                  for (int d=0;d<8;d++) if ((int)short_buffer[d]==local_opt)
                    for (int pos_j=shuffle[d],j=d;j<length && pos_j<m;j+=8,pos_j++)
                      {
                        int t=h0[j];
                        if (t>=opt) 
                          if (t>opt) 
                            opt=t,n_best=1,query_end=i,target_end=min_last_j=max_last_j=pos_j;
                          else 
                            {
                              n_best++;
                              if (direction==1) query_end=i,target_end=pos_j;
                              if (pos_j>max_last_j) max_last_j=pos_j;
                              if (pos_j<min_last_j) min_last_j=pos_j;
                            }
                      }
                }
            }
        }
    }
  if (task_id==3 || task_id==4)
    {
#ifdef LOCAL
      Timer t75(75);
#endif
      query_end=n-1;
      for (int j=0;j<length;j++) if (shuffle[j]<m)
        {
          int t=h[n&1][j];
          if (t>=opt) 
            if (t>opt) 
              opt=t,n_best=1,min_last_j=max_last_j=shuffle[j];
            else 
              {
                n_best++;
                if (shuffle[j]>max_last_j) max_last_j=shuffle[j];
                if (shuffle[j]<min_last_j) min_last_j=shuffle[j];
              }
        }
      target_end=(direction==1)?max_last_j:min_last_j;
    }
  first_j=min_last_j-n-(n*match_score-opt)/(-gap_extension)-2;
  if (first_j<0) first_j=0;
  last_j=max_last_j;
  //sprintf(result,"%d %d %d %d",opt,query_end,target_end,n_best);
  gap_open-=gap_extension;
}

void solve_task1()
{
#ifdef LOCAL
  Counter c1(1);
  Timer t98(98);
#endif
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
    {
#ifdef LOCAL
      Timer t13(13);
#endif
      ca[n]=0;
      uchar s=0;
      for (int i=n-1;i>=0;i--) { a[i]=encode(query[i]); s=(s<<2)|a[i]; ca[i]=s; }
      cb[m]=0;
      memset(vfirst,255,sizeof(vfirst));
      s=0;
      for (int i=m-1;i>=0;i--) 
        { 
          b[i]=encode(target[i]); s=(s<<2)|b[i]; cb[i]=s; 
          if (i<=m-8) vnext[i]=vfirst[s],vfirst[s]=i;
        }
    }
#ifdef LOCAL
  Timer t14(14);
#endif
#ifdef DEBUG
  int cc=0;
#endif
  //int opt=-INF,query_end,target_end,n_best;
  opt=-INF;
  match_score<<=16;
  mismatch_score<<=16;
  gap_open<<=16;
  gap_extension<<=16;
  int gap_open_and_extension=gap_open+gap_extension;
  dpstate *dp[3];
  dp[0]=dp_buffer;
  dp[1]=dp[0]+m+5;
  dp[2]=dp[1]+m+5;
  dp[0][0].idx=dp[1][0].idx=-1;
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate *q0=dp[i&1],*q1=dp[(i+1)&1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      if (opt>est_opt) est_opt=opt;
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<1) baseline=1;
      if (i+7<n && 1>=baseline)
        {
#ifdef LOCAL            
          //Timer t11(11);
#endif
          int pos=1,key=ca[i];
          dpstate *q2=dp[2];
          q0[size0+1].idx=m+2;
          int key_exp=ca[i+4];
          for (int j=vfirst[key];j>=0;j=vnext[j])
            {
#ifdef DEBUG
              ASSERT(j+7<m && cb[j]==key);
#endif
              if (valid0[key_exp^cb[j+4]] && (i==0 || j==0 || a[i-1]!=b[j-1]))
                {                         
                  for (;q0[pos].idx<j;pos++) { q2++; *q2=q0[pos]; }
                  if (q0[pos].idx==j)
                    {
                      q2++; *q2=q0[pos]; pos++;
                    }
                  else
                    {
                      q2++; q2->idx=j; q2->f=j; q2->h=-1; 
                    }
                }
            }
          if (q2-dp[2]>pos-1)
            {
              for (;pos<=size0;pos++) { q2++; *q2=q0[pos]; }
              q0=dp[2];
              size0=q2-dp[2];
            }
        }
      if (size0==0)
        ;
      else if (i==n)
        {
          for (int cl=1;cl<=size0;cl++)
            {
#ifdef DEBUG
              cc++;
#endif
              int j=q0[cl].idx;
              int f0=q0[cl].f;
              int real_score=(f0>>16);
              if (real_score>=opt)
                if (real_score>opt)
                  {
                    opt=real_score; n_best=1; query_end=i; target_end=j;
                    last_j=j-1; first_j=f0&((1<<16)-1);
                  }
                else 
                  {
                    n_best++;
                    if (direction==1) { query_end=i; target_end=j; }
                    if (j-1>last_j) last_j=j-1;
                    int pos_j=f0&((1<<16)-1);
                    if (pos_j<first_j) first_j=pos_j;
                  }
            }
        }
      else
        {
#ifdef LOCAL
          //Timer t12(12);
#endif
#ifdef DEBUG
          for (int cl=1;cl+1<=size0;cl++) ASSERT(q0[cl].idx<q0[cl+1].idx);
          ASSERT(q1==dp[(i+1)&1] && q1->idx==-1);
#endif
          baseline<<=16;
          int g0=-1,ckey=a[i];
          q0[size0+1].idx=m;
          q0[size0+1].f=-1;
          q0++;
          while (q0->idx<m)
            {
#ifdef DEBUG
              cc++;
#endif
              int j=q0->idx;
              int f0=q0->f;
              int h0=q0->h;
              q0++;
              if (g0>f0) f0=g0;
              if (h0>=f0) 
                {
                  int t=h0+gap_extension;
                  if (t>=baseline)
                    {
                      if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                      else if (t>q1->h) q1->h=t;
                    }
                  if (ckey==b[j])
                    {
#ifdef DEBUG
                      ASSERT(h0+match_score>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=h0+match_score; q1->h=-1;
                    }
                  else
                    j++;
                  if (h0+gap_open>=g0) { g0=-1; continue; }
                  g0+=gap_extension;
                }
              else
                {
#ifdef DEBUG
                  ASSERT(f0>=0);
#endif
                  int real_score=(f0>>16);
                  if (real_score>=opt)
                    if (real_score>opt)
                      {
                        opt=real_score; n_best=1; query_end=i; target_end=j;
                        last_j=j-1; first_j=f0&((1<<16)-1);
                      }
                    else 
                      {
                        n_best++;
                        if (direction==1) { query_end=i; target_end=j; }
                        if (j-1>last_j) last_j=j-1;
                        int pos_j=f0&((1<<16)-1);
                        if (pos_j<first_j) first_j=pos_j;
                      }
                  g0+=gap_extension;
                  if (ckey==b[j])
                    {
                      int t=h0+gap_extension;
                      if (t>=baseline)
                        {
                          if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                          else if (t>q1->h) q1->h=t;
                        }
#ifdef DEBUG
                      ASSERT(f0+match_score>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=f0+match_score; q1->h=-1;
                    }
                  else
                    {
                      int t=f0+gap_open_and_extension;
                      if (t>=baseline)
                        {
                          if (t>g0) g0=t;
                          if (h0+gap_extension>t) t=h0+gap_extension;
                          if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                          else if (t>q1->h) q1->h=t;
                        }
                      else
                        {
                          t=h0+gap_extension;
                          if (t>=baseline)
                            if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                            else if (t>q1->h) q1->h=t;
                        }
                      t=f0+mismatch_score;
                      j++;
                      if (t>=baseline) { q1++; q1->idx=j; q1->f=t; q1->h=-1; }
                    }
                }
              if (g0<baseline) { g0=-1; continue; }
              int nextj=q0->idx;
              while (j<nextj)
                {
#ifdef DEBUG
                  cc++;
                  ASSERT(g0>0);
#endif
                  if (ckey==b[j])
                    {
                      int t=g0+match_score;
#ifdef DEBUG
                      ASSERT(t>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=t; q1->h=-1;
                    }
                  else
                    j++;
                  g0+=gap_extension;
                  if (g0<baseline) { g0=-1; break; }
                }
            }
          if (q0->f>=0)
            {
              int j=q0->idx;
              int f0=q0->f;
              int real_score=(f0>>16);
              if (real_score>=opt)
                if (real_score>opt)
                  {
                    opt=real_score; n_best=1; query_end=i; target_end=j;
                    last_j=j-1; first_j=f0&((1<<16)-1);
                  }
                else 
                  {
                    n_best++;
                    if (direction==1) { query_end=i; target_end=j; }
                    if (j-1>last_j) last_j=j-1;
                    int pos_j=f0&((1<<16)-1);
                    if (pos_j<first_j) first_j=pos_j;
                  }
            }
          size0=q1-dp[(i+1)&1];
        }
    }
  match_score>>=16;
  mismatch_score>>=16;
  gap_open>>=16;
  gap_extension>>=16;
#ifdef DEBUG
  printf("ratio = %.3lf    opt = %d\n",(double)cc/(double)(n*m),opt);
#endif
  //sprintf(result,"%d %d %d %d",opt,query_end-1,target_end-1,n_best);
  query_end--; target_end--;
}

void solve_task2()
{
#ifdef LOCAL
  Counter c2(2);
#endif
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
  int gap_open_and_extension=gap_open+gap_extension;
    {
#ifdef LOCAL
      Timer t23(23);
#endif
      if (counter_task2==60000) counter_task2=0;
      if (counter_task2==0) memset(dp_large,0,sizeof(dp_large));
      counter_task2++;
      uchar s=0;
      ca[n]=0;
      for (int i=n-1;i>=0;i--) { a[i]=encode(query[n-1-i]); s=(s<<2)|a[i]; ca[i]=s; }
      cb[m]=0;
      s=0;
      for (int i=m-1;i>=0;i--) { b[i]=encode(target[m-1-i]); s=(s<<2)|b[i]; cb[i]=s; }
      memset(vfirst,255,64*sizeof(int));
      for (int i=m-7;i>=0;i--)
        {
          int set=cb[i]&63;
          vnext[i]=vfirst[set];
          vfirst[set]=i;
        }
      if (n>=8) for (int src=0;src+7<m;src++) if (ca[0]==cb[src] && ca[4]==cb[src+4])
        {
          bool escape=false;
          int i=0,j=src,cost=0,greedy_cost=0;
          int e2=gap_open_and_extension+gap_extension;
          int e3=e2+gap_extension;
          int e4=e3+gap_extension;
          for (;i<n;)
            {
              if (cost<0) { escape=true; cost=0; }
              if (j>=m) { cost+=gap_open-((n-i)<<1)/**gap_extension*/; break; }
              else if (a[i]==b[j]) { cost+=match_score; i++; j++; }
              else
                {
                  int t=cost+gap_open-((n-i)<<1);//*gap_extension;
                  if (t>greedy_cost) greedy_cost=t;
                  if (ca[i]==cb[j+1]) { cost+=gap_open_and_extension; j++; continue; }
                  if (ca[i+1]==cb[j]) { cost+=gap_open_and_extension; i++; continue; }
                  if (ca[i]==cb[j+2]) { cost+=e2; j+=2; continue; }
                  if (ca[i+2]==cb[j]) { cost+=e2; i+=2; continue; }
                  if (ca[i]==cb[j+3]) { cost+=e3; j+=3; continue; }
                  if (ca[i+3]==cb[j]) { cost+=e3; i+=3; continue; }
                  if (ca[i]==cb[j+4]) { cost+=e4; j+=4; continue; }
                  if (ca[i+4]==cb[j]) { cost+=e4; i+=4; continue; }
                  cost+=mismatch_score; i++; j++;
                }
            }
          if (cost>greedy_cost) greedy_cost=cost;
          if (greedy_cost>est_opt) est_opt=greedy_cost;
          if (greedy_cost==n)//*match_score)
            {
              first_j=src;
              last_j=src+n-1;
              int n_best=1,target_end=first_j;
              for (int j=src+1;j+n-1<m;j++) if (ca[0]==cb[j] && memcmp(a,b+j,n*sizeof(int))==0)
                {
                  last_j=j+n-1;
                  n_best++;
                  if (direction==0) target_end=j;
                }
              target_end=m-1-target_end;
              swap(first_j,last_j);
              first_j=m-1-first_j;
              last_j=m-1-last_j;
              //sprintf(result,"%d %d %d %d",greedy_cost,n-1,last_j,n_best);
              opt = greedy_cost;
              query_end=n-1;
              target_end=last_j;
              return;
            }
        }
    }
  match_score<<=16;
  mismatch_score<<=16;
  gap_open<<=16;
  gap_extension<<=16;
  gap_open_and_extension<<=16;
#ifdef DEBUG
  int cc=0;
#endif
  //int opt=-INF,query_end,target_end,n_best;
  opt=-INF;
  int best_end_points[MAX_LENGTH];
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate_large *dp0=dp_large[i],*dp1=dp_large[i+1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      int *q0=q[i&1],*q1=q[(i+1)&1];
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<0) baseline=0;
      if (i+6<n)
        {
#ifdef LOCAL            
          Timer t21(21);
#endif
          int pos=0,*q2=q[2],key_exp=ca[i+3];
          q0[size0]=m+2;
          for (int j=vfirst[ca[i]&((1<<6)-1)];j>=0;j=vnext[j])
            {
#ifdef DEBUG
              ASSERT(j+6<m && (cb[j]&((1<<6)-1))==ca[i]&((1<<6)-1));
#endif
              if (valid0[key_exp^cb[j+3]] && (i==0 || j==0 || a[i-1]!=b[j-1]))
                {
                  for (;q0[pos]<j;pos++) { *q2=q0[pos]; q2++; }
                  if (q0[pos]==j)
                    {
                      *q2=q0[pos]; q2++; pos++;
                    }
                  else
                    {
                      dp0[j].f=counter_task2; dp0[j].h=-1; *q2=j; q2++;
                    }
                }
            }
          if (q2-q[2]>pos)
            {
              for (;pos<size0;pos++) { *q2=q0[pos]; q2++; }
              q0=q[2];
              size0=q2-q[2];
            }
        }
      if (size0==0)
        ;
      else if (i==n)
        {
          for (int cl=0;cl<size0;cl++)
            {
              int j=q0[cl];
              int f0=dp0[j].f;
              int h0=dp0[j].h;
              if (h0>f0) f0=h0;
              dp0[j].f=f0; 
              int real_score=(f0>>16);
              if (real_score>=opt)
                if (real_score>opt)
                  {
                    opt=real_score;
                    n_best=0;
                    best_end_points[n_best++]=j;
                  }
                else
                  best_end_points[n_best++]=j;
            }
        }
      else
        {
#ifdef LOCAL
          Timer t22(22);
#endif
          int g0=-1,last_added=-1;
          int ckey=a[i];
          q0[size0]=m;
          baseline<<=16;
          int reduce=-5-((n-i)<<1);//((gap_open+(n-i)*gap_extension)>>16);
          int local_est_opt=0;
          for (int cl=0;cl<size0;cl++)
            {
#ifdef DEBUG
              cc++;
#endif
              int j=q0[cl];
              int f0=dp0[j].f;
              int h0=dp0[j].h;
              if (g0>f0) f0=g0;
              g0+=gap_extension;
              if (h0>=f0) f0=h0;
              dp0[j].f=f0;
#ifdef DEBUG
              ASSERT(f0>=0);
#endif
              int real_score=(f0>>16);
              if (real_score>local_est_opt) local_est_opt=real_score;
              if (j<m)
                if (ckey==b[j])
                  {
                    int t=h0+gap_extension;
                    if (t>=baseline)
                      {
                        if (j!=last_added) { last_added=j; dp1[j].f=-1; dp1[j].h=t; *q1=j; q1++; }
                        else if (t>dp1[j].h) dp1[j].h=t;
                      }
#ifdef DEBUG
                    ASSERT(f0+match_score>=0);
#endif
                    j++; last_added=j; dp1[j].f=f0+match_score; dp1[j].h=-1; *q1=j; q1++;
                  }
                else
                  {
                    int t=f0+gap_open_and_extension;
                    if (t>=baseline)
                      {
                        if (t>g0) g0=t;
                        if (h0+gap_extension>t) t=h0+gap_extension;
                        if (j!=last_added) { last_added=j; dp1[j].f=-1; dp1[j].h=t; *q1=j; q1++; }
                        else if (t>dp1[j].h) dp1[j].h=t;
                      }
                    else
                      {
                        t=h0+gap_extension;
                        if (t>=baseline)
                          if (j!=last_added) { last_added=j; dp1[j].f=-1; dp1[j].h=t; *q1=j; q1++; }
                          else if (t>dp1[j].h) dp1[j].h=t;
                      }
                    t=f0+mismatch_score;
                    j++;
                    if (t>=baseline) { last_added=j; dp1[j].f=t; dp1[j].h=-1; *q1=j; q1++; }
                  }
              else
                {
                  int t=max(f0+gap_open_and_extension,h0+gap_extension);
                  if (t>=baseline)
                    {
                      if (t>g0) g0=t;
                      if (j!=last_added) { last_added=j; dp1[j].f=-1; dp1[j].h=t; *q1=j; q1++; }
                      else if (t>dp1[j].h) dp1[j].h=t;
                    }
                }
              if (g0<baseline) { g0=-1; continue; }
              int nextj=q0[cl+1];
              while (j<nextj)
                {
#ifdef DEBUG
                  cc++;
                  ASSERT(g0>=0);
#endif
                  dp0[j].f=g0;
                  dp0[j].h=-1;
                  if (ckey==b[j])
                    {
                      int t=g0+match_score;
#ifdef DEBUG
                      ASSERT(t>=baseline);
#endif
                      j++; last_added=j; dp1[j].f=t; dp1[j].h=-1; *q1=j; q1++;
                    }
                  else
                    j++;
                  g0+=gap_extension;
                  if (g0<baseline) { g0=-1; break; }
                }
            }
          local_est_opt+=reduce;
          if (local_est_opt>est_opt) est_opt=local_est_opt;
          size0=q1-q[(i+1)&1];
        }
    }
    {
#ifdef LOCAL
      Timer t24(24);
#endif
#ifdef DEBUG
      for (int i=0;i+1<n_best;i++) ASSERT(best_end_points[i]<best_end_points[i+1]);
#endif
      first_j=m-best_end_points[n_best-1];
      reverse(best_end_points,best_end_points+n_best);
      for (int i=0;i<n_best;i++)
        { 
          int pos=best_end_points[i];
          q[n&1][i]=pos;
          vf[n][pos]=true;
          vh[n][pos]=false; 
        }
      sizeq[n&1]=n_best;
      vector<ipair> respos_all;
      last_j=0;
      for (int i=n;i>=0;i--)
        {
          int size0=sizeq[i&1],last_added=-1,last_g=INF;
          int *q0=q[i&1],*q1=q[(i-1)&1];
          bool *vh0=vh[i],*vf0=vf[i],*vh1=vh[i-1],*vf1=vf[i-1];
#ifdef DEBUG
          for (int cl=0;cl+1<size0;cl++) ASSERT(q0[cl]>q0[cl+1]);
#endif
          for (int cl=0;cl<size0;cl++)
            {
              int j=q0[cl];
              if (vf0[j] && dp_large[i][j].f==dp_large[i][j].h) vh0[j]=true;
              if (vh0[j])
                {
#ifdef DEBUG
                  ASSERT(dp_large[i][j].h>=0);
#endif                
                  if (i>0 && dp_large[i-1][j].h+gap_extension==dp_large[i][j].h)
                    {
                      if (last_added!=j) { vf1[j]=false; last_added=j; *q1=j; q1++; }
                      vh1[j]=true;
                    }
                  if (i>0 && dp_large[i-1][j].f+gap_open_and_extension==dp_large[i][j].h)
                    {
                      if (last_added!=j) { vh1[j]=false; last_added=j; *q1=j; q1++; }
                      vf1[j]=true;
                    }
                }
              if (vf[i][j])
                {
#ifdef DEBUG
                  ASSERT(dp_large[i][j].f>=0);
#endif
                  if ((dp_large[i][j].f>>16)==0)
                    {
                      respos_all.push_back(MP(n-i,m-j));
                      int pos_j=m-1-j;
                      if (pos_j>last_j) last_j=pos_j; 
                    }
                  if (i>0)
                    {
                      if (dp_large[i-1][j].f+gap_open_and_extension==dp_large[i][j].f)
                        {
                          if (last_added!=j) { vh1[j]=false; last_added=j; *q1=j; q1++; }
                          vf1[j]=true;
                        }
                      if (dp_large[i-1][j].h+gap_extension==dp_large[i][j].f)
                        {
                          if (last_added!=j) { vf1[j]=false; last_added=j; *q1=j; q1++; }
                          vh1[j]=true;
                        }
                    }
                  if (i>0 && j>0 && 
                      dp_large[i-1][j-1].f+((a[i-1]==b[j-1])?match_score:mismatch_score)==dp_large[i][j].f)
                    {
                      if (last_added!=j-1) { vh1[j-1]=false; last_added=j-1; *q1=j-1; q1++; }
                      vf1[j-1]=true;
                    }
                  int t=dp_large[i][j].f+gap_extension;
                  if (t<last_g) last_g=t;
                }
              int nextj=(cl+1<size0)?q0[cl+1]:-1;
              while (1)
                {
                  last_g-=gap_extension;
                  j--;
                  if (j==nextj) break;
                  if ((dp_large[i][j].f&((1<<16)-1))!=counter_task2
                      && (dp_large[i][j].h&((1<<16)-1))!=counter_task2) 
                    { 
                      last_g=INF; break; 
                    }
                  if (dp_large[i][j].f+gap_open_and_extension==last_g)
                    {
                      vh0[j]=false;
                      vf0[j]=true;
                      q0[cl--]=j;
                      break;
                    }
                }
            }
          sizeq[(i-1)&1]=q1-q[(i-1)&1];
        }
      n_best=SIZE(respos_all);
#ifdef DEBUG
      ASSERT(n_best>0);
#endif
      if (direction==0)
        {
          query_end=respos_all.begin()->first;
          target_end=respos_all.begin()->second;
        }
      else
        {
          query_end=(--respos_all.end())->first;
          target_end=(--respos_all.end())->second;
        }
#ifdef DEBUG
      printf("ratio = %.3lf    opt = %d\n",(double)cc/(double)(n*m),opt);
#endif
      //sprintf(result,"%d %d %d %d",opt,query_end-1,target_end-1,n_best);
      query_end--;
      target_end--;
    }
  match_score>>=16;
  mismatch_score>>=16;
  gap_open>>=16;
  gap_extension>>=16;
}

void solve_task3()
{
#ifdef LOCAL
  Counter c3(3);
  Timer t98(98);
#endif
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
  int gap_open_and_extension=gap_open+gap_extension;
    {
#ifdef LOCAL
      Timer t33(33);
#endif
      ca[n]=0;
      uchar s=0;
      for (int i=n-1;i>=0;i--) { a[i]=encode(query[i]); s=(s<<2)|a[i]; ca[i]=s; }
      cb[m]=0;
      memset(vfirst,255,sizeof(vfirst));
      s=0;
      for (int i=m-1;i>=0;i--) 
        { 
          b[i]=encode(target[i]); s=(s<<2)|b[i]; cb[i]=s; 
          if (i<=m-8) vnext[i]=vfirst[s],vfirst[s]=i;
        }
      if (n>=8) for (int src=0;src+7<m;src++) if (ca[0]==cb[src] && ca[4]==cb[src+4])
        {
          int i=0,j=src,cost=0,greedy_cost=0;
          int e2=gap_open_and_extension+gap_extension;
          int e3=e2+gap_extension;
          int e4=e3+gap_extension;
          bool escape=false;
          for (;i<n;)
            {
              if (cost<0) { escape=true; cost=0; }
              if (j>=m) { cost+=gap_open-((n-i)<<1)/**gap_extension*/; break; }
              else if (a[i]==b[j]) { cost+=match_score; i++; j++; }
              else
                {
                  int t=cost+gap_open-((n-i)<<1);//*gap_extension;
                  if (t>greedy_cost) greedy_cost=t;
                  if (ca[i]==cb[j+1]) { cost+=gap_open_and_extension; j++; continue; }
                  if (ca[i+1]==cb[j]) { cost+=gap_open_and_extension; i++; continue; }
                  if (ca[i]==cb[j+2]) { cost+=e2; j+=2; continue; }
                  if (ca[i+2]==cb[j]) { cost+=e2; i+=2; continue; }
                  if (ca[i]==cb[j+3]) { cost+=e3; j+=3; continue; }
                  if (ca[i+3]==cb[j]) { cost+=e3; i+=3; continue; }
                  if (ca[i]==cb[j+4]) { cost+=e4; j+=4; continue; }
                  if (ca[i+4]==cb[j]) { cost+=e4; i+=4; continue; }
                  cost+=mismatch_score; i++; j++;
                }
            }
          if (cost>greedy_cost) greedy_cost=cost;
          if (greedy_cost>est_opt) est_opt=greedy_cost;
          if (greedy_cost==n)//*match_score)
            {
              first_j=src;
              last_j=src+n-1;
              int n_best=1,target_end=last_j;
              for (int j=src+1;j+n-1<m;j++) if (ca[0]==cb[j] && memcmp(query,target+j,n*sizeof(char))==0)
                {
                  last_j=j+n-1;
                  n_best++;
                  if (direction==1) target_end=last_j;
                }
              //printf(result,"%d %d %d %d",greedy_cost,n-1,target_end,n_best);
              opt = greedy_cost;
              query_end = n-1;
              return;
            }
        }
    }
#ifdef DEBUG
  int cc=0;
#endif
  int opt=-INF,target_end,n_best;
  opt=-INF;
  match_score<<=16;
  mismatch_score<<=16;
  gap_open<<=16;
  gap_extension<<=16;
  gap_open_and_extension<<=16;
  dpstate *dp[3];
  dp[0]=dp_buffer;
  dp[1]=dp[0]+m+5;
  dp[2]=dp[1]+m+5;
  dp[0][0].idx=dp[1][0].idx=-1;
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate *q0=dp[i&1],*q1=dp[(i+1)&1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<1) baseline=1;
      if (i+7<n && 1>=baseline)
        {
#ifdef LOCAL            
          Timer t31(31);
#endif
          int pos=1,key=ca[i];
          dpstate *q2=dp[2];
          q0[size0+1].idx=m+2;
          int key_exp=ca[i+4];
          for (int j=vfirst[key];j>=0;j=vnext[j])
            {
#ifdef DEBUG
              ASSERT(j+7<m && cb[j]==key);
#endif
              if (valid0[key_exp^cb[j+4]] && (i==0 || j==0 || a[i-1]!=b[j-1]))
                {                         
                  for (;q0[pos].idx<j;pos++) { q2++; *q2=q0[pos]; }
                  if (q0[pos].idx==j)
                    {
                      q2++; *q2=q0[pos]; pos++;
                    }
                  else
                    {
                      q2++; q2->idx=j; q2->f=j; q2->h=-1; 
                    }
                }
            }
          if (q2-dp[2]>pos-1)
            {
              for (;pos<=size0;pos++) { q2++; *q2=q0[pos]; }
              q0=dp[2];
              size0=q2-dp[2];
            }
        }
      if (size0==0)
        ;
      else if (i==n)
        {
          for (int cl=1;cl<=size0;cl++)
            {
              int j=q0[cl].idx;
              int f0=max(q0[cl].f,q0[cl].h);
              int real_score=(f0>>16);
              if (real_score>=opt)
                if (real_score>opt)
                  {
                    opt=real_score; n_best=1; target_end=j;
                    last_j=j-1; first_j=f0&((1<<16)-1);
                  }
                else
                  {
                    if (j-1!=last_j) n_best++;
                    if (direction==0) { if (j<target_end) target_end=j; }
                    else { if (j>target_end) target_end=j; }
                    if (j-1>last_j) last_j=j-1;
                    int pos_j=f0&((1<<16)-1); if (pos_j<first_j) first_j=pos_j;
                  }
            }
        }
      else
        {
#ifdef LOCAL
          Timer t32(32);
#endif
#ifdef DEBUG
          for (int cl=1;cl+1<=size0;cl++) ASSERT(q0[cl].idx<q0[cl+1].idx);
          ASSERT(q1==dp[(i+1)&1] && q1->idx==-1);
#endif
          baseline<<=16;
          int g0=-1,ckey=a[i];
          int reduce=-5-((n-i)<<1);//((gap_open+(n-i)*gap_extension)>>16);
          int local_est_opt=0;
          if (q0[size0].idx==m)
            {
              int f0=max(q0[size0].f,q0[size0].h-gap_open);
              int real_score=(f0>>16)+reduce;
              if (real_score>=0)
                {
                  int j=m;
                  if (real_score>=opt) 
                    if (real_score>opt)
                      {
                        opt=real_score; n_best=1; target_end=j;
                        last_j=j-1; first_j=f0&((1<<16)-1);
                      }
                    else 
                      {
                        if (j-1!=last_j) n_best++;
                        if (direction==0) { if (j<target_end) target_end=j; }
                        else { if (j>target_end) target_end=j; }
                        if (j-1>last_j) last_j=j-1;
                        int pos_j=f0&((1<<16)-1);
                        if (pos_j<first_j) first_j=pos_j;
                      }
                  if (real_score>est_opt) est_opt=real_score;
                }
            }
          else
            q0[size0+1].idx=m;
          q0++;
          while (q0->idx<m)
            {
#ifdef DEBUG
              cc++;
#endif
              int j=q0->idx;
              int f0=q0->f;
              int h0=q0->h;
              q0++;
              if (g0>f0) f0=g0;
              if (h0>=f0) 
                {
                  int t=h0+gap_extension;
                  if (t>=baseline)
                    {
                      if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                      else if (t>q1->h) q1->h=t;
                    }
                  if (ckey==b[j])
                    {
#ifdef DEBUG
                      ASSERT(h0+match_score>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=h0+match_score; q1->h=-1;
                    }
                  else
                    j++;
                  if (h0+gap_open>=g0) { g0=-1; continue; }
                  g0+=gap_extension;
                }
              else
                {
#ifdef DEBUG
                  ASSERT(f0>=0);
#endif
                  int real_score=(f0>>16);
                  if (real_score>local_est_opt) local_est_opt=real_score;
                  g0+=gap_extension;
                  if (ckey==b[j])
                    {
                      int t=h0+gap_extension;
                      if (t>=baseline)
                        {
                          if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                          else if (t>q1->h) q1->h=t;
                        }
#ifdef DEBUG
                      ASSERT(f0+match_score>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=f0+match_score; q1->h=-1;
                    }
                  else
                    {
                      int t=f0+gap_open_and_extension;
                      if (t>=baseline)
                        {
                          if (t>g0) g0=t;
                          if (h0+gap_extension>t) t=h0+gap_extension;
                          if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                          else if (t>q1->h) q1->h=t;
                        }
                      else
                        {
                          t=h0+gap_extension;
                          if (t>=baseline)
                            if (j!=q1->idx) { q1++; q1->idx=j; q1->f=-1; q1->h=t; }
                            else if (t>q1->h) q1->h=t;
                        }
                      t=f0+mismatch_score;
                      j++;
                      if (t>=baseline) { q1++; q1->idx=j; q1->f=t; q1->h=-1; }
                    }
                }
              if (g0<baseline) { g0=-1; continue; }
              int nextj=q0->idx;
              while (j<nextj)
                {
#ifdef DEBUG
                  cc++;
                  ASSERT(g0>0);
#endif
                  if (ckey==b[j])
                    {
                      int t=g0+match_score;
#ifdef DEBUG
                      ASSERT(t>0);
#endif
                      j++; q1++; q1->idx=j; q1->f=t; q1->h=-1;
                    }
                  else
                    j++;
                  g0+=gap_extension;
                  if (g0<baseline) { g0=-1; break; }
                }
            }
          local_est_opt+=reduce;
          if (local_est_opt>est_opt) est_opt=local_est_opt;
          size0=q1-dp[(i+1)&1];
        }
    }
  match_score>>=16;
  mismatch_score>>=16;
  gap_open>>=16;
  gap_extension>>=16;
#ifdef DEBUG
  printf("ratio = %.3lf    opt = %d\n",(double)cc/(double)(n*m),opt);
#endif
  //sprintf(result,"%d %d %d %d",opt,n-1,target_end-1,n_best);
  query_end=n-1;
  target_end--;
}

void solve()
{
  if (measure_id!=0) { brute_force(); return; }
  if (task_id==4) { brute_force(); return; }
  if (task_id==1)
    {
      solve_task1();
      return;
    }
  if (task_id==2)
    {
      solve_task2();
      return;
    }
  if (task_id==3)
    {
      solve_task3();
      return;
    }
}

#ifdef LOCAL

char __buffer[128<<10];
int __bufferL=0;
int __current=0;
FILE *fin;

char __readchar()
{
  if (__current==__bufferL)
    {
      __bufferL=fread(__buffer,sizeof(char),128<<10,fin);
      if (__bufferL==0) return 0;
      __current=0;
    }
  return __buffer[__current++];
}

void __readstring(char *s)
{
  char c;
  do{ c=__readchar(); }while (c<=32);
  for (;c>32;c=__readchar()) { *s=c; s++; }
  *s=0;
}

char __buffer_used_for_nextint[1024];

int __nextint()
{
  __readstring(__buffer_used_for_nextint);
  int op=1;
  char *p=__buffer_used_for_nextint;
  if (*p=='-') p++,op=-1;
  int val=0;
  for (;*p!=0;p++) val=val*10+(*p-'0');
  if (op<0) val=-val;
  return val;
}

int main(int argc,char **args)
{
  string filename="..\\Data\\1.txt";
  if (argc==2) filename=args[1];
  fin=fopen(filename.c_str(),"r");
  uint64 start_time=rdtsc();
  int case_id;
  int c[2000];
  memset(c,0,sizeof(c));
  for (case_id=0;case_id<50000;case_id++)
    {
      char target[1500],query[1500];
      __readstring(target);
      if (target[0]=='#') break;
      __readstring(query);
      int query_start_clip=__nextint();
      int query_end_clip=__nextint();
      int match_score=__nextint();
      int mismatch_score=__nextint();
      int gap_open=__nextint();
      int gap_extension=__nextint();
      int direction=__nextint();
      int opt=__nextint();
      int query_end=__nextint();
      int target_end=__nextint();
      int n_best=__nextint();
      //if (case_id!=29302) continue;
#ifdef DEBUG
      printf("Case : %d\n",case_id);
#endif
      __opt=opt;
      c[opt/20]++;
      string ret=process(target,query,query_start_clip,query_end_clip,match_score,mismatch_score,gap_open,gap_extension,direction);
      int r_opt,r_query_end,r_target_end,r_n_best;
      sscanf(ret.c_str(),"%d%d%d%d",&r_opt,&r_query_end,&r_target_end,&r_n_best);
      if (r_opt!=opt || r_query_end!=query_end || r_target_end!=target_end || r_n_best!=n_best)
        {
          printf("ERROR in CASE # %d\n",case_id);
          return 0;
        }
    }
  uint64 end_time=rdtsc();
  fclose(fin);
  printf("ALLPASSED   %d  : %.3lf\n",case_id,(double)(end_time-start_time)/1e9);
  Counter::show();
  Timer::show();
  //for (int i=0;i<2000;i++) if (c[i]>0) printf("%d %d\n",i,c[i]);
  //printf("n_states = %d\n",n_states);
  printf("cc1 = %d\n",cc1);
  printf("cc2 = %d\n",cc2);
  printf("cc3 = %d\n",cc3);
  return 0;
}

#endif
