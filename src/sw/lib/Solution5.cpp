// Coder: ACRush
// Submission: 27
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=19849563&subnum=27
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

#include <emmintrin.h>
#include <pmmintrin.h>
#include <xmmintrin.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "Solution5.h"

using namespace std;

#define SIZE(A) ((int)A.size())
#define LENGTH(A) ((int)A.length())
typedef int64_t int64;
typedef uint64_t uint64;
typedef uint32_t uint;
typedef uint16_t ushort;
typedef uint8_t uchar;
typedef pair<int,int> ipair;
#define MP(A,B) make_pair(A,B)
typedef vector<string> VS;
typedef vector<int> VI;
#define two(X) (1<<(X))
#define contain(S,X) ((S&two(X))>0)
const double pi=acos(-1.0);
inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }

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

void Solution5::init()
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

int get_task_id(int qsc,int qec)
{
  if (qsc==1 && qec==1) return 1;
  if (qsc==0 && qec==1) return 2;
  if (qsc==1 && qec==0) return 3;
  if (qsc==0 && qec==0) return 4;
  return 0;
}

void Solution5::resize(int mem) {
    int i, prev;
    if(mem < MAX_LENGTH) return;
    prev = MAX_LENGTH;
    MAX_LENGTH = mem;
    tmap_roundup32(MAX_LENGTH);
    target = (char*)tmap_realloc(target, sizeof(char) * MAX_LENGTH, "target");
    query = (char*)tmap_realloc(query, sizeof(char) * MAX_LENGTH, "query");
    a = (int*)tmap_realloc(a, sizeof(int) * MAX_LENGTH, "a");
    b = (int*)tmap_realloc(b, sizeof(int) * MAX_LENGTH, "b");
    ca = (int*)tmap_realloc(ca, sizeof(int) * MAX_LENGTH, "ca");
    cb = (int*)tmap_realloc(cb, sizeof(int) * MAX_LENGTH, "cb");
    dp_buffer = (dpstate_t*)tmap_realloc(dp_buffer, sizeof(dpstate_t) * (100 + (MAX_LENGTH * 3)), "dp_buffer");
    for(i=0;i<3;i++) q[i] = (int*)tmap_realloc(q[i], sizeof(int) * MAX_LENGTH, "q[i]");
    dp_large = (dpstate_large_t**)tmap_realloc(dp_large, sizeof(dpstate_large_t*) * MAX_LENGTH, "dp_large");
    vf = (bool**)tmap_realloc(vf, sizeof(bool*) * MAX_LENGTH, "vf");
    vh = (bool**)tmap_realloc(vh, sizeof(bool*) * MAX_LENGTH, "vh");
    for(i=0;i<MAX_LENGTH;i++) {
        if(prev <= i) {
            dp_large[i] = NULL;
            vf[i] = vh[i] = NULL;
        }
        dp_large[i] = (dpstate_large_t*)tmap_realloc(dp_large[i], sizeof(dpstate_large_t) * MAX_LENGTH, "dp_large");
        vf[i] = (bool*)tmap_realloc(vf[i], sizeof(bool) * MAX_LENGTH, "vf[i]");
        vh[i] = (bool*)tmap_realloc(vh[i], sizeof(bool) * MAX_LENGTH, "vh");
    }
    hbuffer16 = (short*)tmap_realloc(hbuffer16, sizeof(short) * (100 + (MAX_LENGTH * 3)), "hbuffer16");
    hbuffer8 = (uint8_t*)tmap_realloc(hbuffer8, sizeof(uint8_t) * (100 + (MAX_LENGTH * 3)), "hbuffer8");
    for(i=0;i<256;i++) vgraph[i] = (int*)tmap_realloc(vgraph[i], sizeof(int) * MAX_LENGTH, "vgraph[i]");
    vnext = (int*)tmap_realloc(vnext, sizeof(int) * MAX_LENGTH, "vnext");
    shuffle = (int*)tmap_realloc(shuffle, sizeof(int) * MAX_LENGTH, "shuffle");
    for(i=0;i<4;i++) {
        cost8[i] = (uint8_t*)tmap_realloc(cost8[i], sizeof(uint8_t) * MAX_LENGTH, "cost8[i]");
        cost16[i] = (short*)tmap_realloc(cost16[i], sizeof(short) * MAX_LENGTH, "cost16");
    }
    transfer_h = (short*)tmap_realloc(transfer_h, sizeof(short) * MAX_LENGTH, "transfer_h");
    transfer_e = (short*)tmap_realloc(transfer_e, sizeof(short) * MAX_LENGTH, "transfer_e");
    best_end_points = (int*)tmap_realloc(best_end_points, sizeof(int) * MAX_LENGTH, "best_end_points");
}

Solution5::Solution5() {
    int i;
    n_best = 1;
    m = n = 0;
    testcase_id=0;
    MAX_LENGTH=0;
    target = query = NULL;
    a = b = ca = cb = NULL;
    dp_buffer= NULL;
    for(i=0;i<3;i++) q[i] = NULL;
    dp_large = NULL;
    vf = vh = NULL;
    hbuffer16 = NULL;
    hbuffer8 = NULL;
    for(i=0;i<256;i++) vgraph[i] = NULL;
    vnext = NULL;
    shuffle = NULL;
    for(i=0;i<4;i++) {
        cost8[i] = NULL;
        cost16[i] = NULL;
    }
    transfer_h = transfer_e = NULL;
    best_end_points = NULL;

    // resize
    resize(1040);
  
    max_qlen = 512;
    max_tlen = 1024;
}

Solution5::~Solution5() {
    int i;
    free(target);
    free(query);
    free(a);
    free(b);
    free(ca);
    free(cb);
    free(dp_buffer);
    for(i=0;i<3;i++) free(q[i]);
    for(i=0;i<MAX_LENGTH;i++) {
        free(dp_large[i]);
        free(vf[i]);
        free(vh[i]);
    }
    free(dp_large);
    free(vf);
    free(vh);
    free(hbuffer16);
    free(hbuffer8);
    for(i=0;i<256;i++) free(vgraph[i]);
    free(vnext);
    free(shuffle);
    for(i=0;i<4;i++) {
        free(cost8[i]);
        free(cost16[i]);
    }
    free(transfer_h);
    free(transfer_e);
    free(best_end_points);
}

int Solution5::
process(const string& s_target, const string& s_query, int qsc, int qec, int v_match_score, 
        int v_mismatch_score, int v_gap_open, int v_gap_extension, int v_direction,
        int *_opt, int *_te, int *_qe, int *_n_best)
{
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
  if(l_query < l_target) resize(l_target);
  else resize(l_query);
  est_opt=0;
  if (testcase_id>0 && 
      direction==v_direction && 
      task_id==v_task_idx &&
      n==l_query && 
      m>=l_target &&
      memcmp(v_query,query,n*sizeof(char))==0)
    {
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
#define encode(c) encode_idx[(int)c]

void Solution5::shift_right16(__m128i &a,short first)
{
  __m128i *p=(__m128i*)short_buffer;
  _mm_store_si128(p,a);
  for (int i=7;i>0;i--) short_buffer[i]=short_buffer[i-1];
  short_buffer[0]=first;
  a=_mm_load_si128(p);
}

void Solution5::shift_right8(__m128i &a,uchar first)
{
  __m128i *p=(__m128i*)uchar_buffer;
  _mm_store_si128(p,a);
  for (int i=15;i>0;i--) uchar_buffer[i]=uchar_buffer[i-1];
  uchar_buffer[0]=first;
  a=_mm_load_si128(p);
}

#define ZERO (45)
#define LBOUND1 (30)
#define LBOUND2 (15)

void Solution5::brute_force()
{
    {
      for (int i=0;i<n;i++) a[i]=encode(query[i]);
      for (int j=0;j<m;j++) b[j]=encode(target[j]);
    }
  gap_open+=gap_extension;
  //int opt=1,min_last_j,max_last_j,n_best,query_end,target_end;
  int min_last_j,max_last_j; 
  opt=1;
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
                  //int pos=shuffle[j];
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

void Solution5::solve_task1()
{
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
    {
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
  dpstate_t *dp[3];
  dp[0]=dp_buffer;
  dp[1]=dp[0]+m+5;
  dp[2]=dp[1]+m+5;
  dp[0][0].idx=dp[1][0].idx=-1;
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate_t *q0=dp[i&1],*q1=dp[(i+1)&1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      if (opt>est_opt) est_opt=opt;
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<1) baseline=1;
      if (i+7<n && 1>=baseline)
        {
          int pos=1,key=ca[i];
          dpstate_t *q2=dp[2];
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

void Solution5::solve_task2()
{
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
  int gap_open_and_extension=gap_open+gap_extension;
    {
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
              n_best=1;
              //int target_end=first_j;
              target_end=first_j;
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
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate_large_t *dp0=dp_large[i],*dp1=dp_large[i+1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      int *q0=q[i&1],*q1=q[(i+1)&1];
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<0) baseline=0;
      if (i+6<n)
        {
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
#ifdef DEBUG
      for (int i=0;i+1<n_best;i++) ASSERT(best_end_points[i]<best_end_points[i+1]);
#endif
      // FIXME: (n_best < 1) can occur
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

void Solution5::solve_task3()
{
#ifdef DEBUG
  ASSERT(measure_id==0);
#endif
  int gap_open_and_extension=gap_open+gap_extension;
    {
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
              n_best=1;
              //int target_end=last_j;
              target_end=last_j;
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
  //int opt=-INF,target_end,n_best;
  opt=-INF;
  match_score<<=16;
  mismatch_score<<=16;
  gap_open<<=16;
  gap_extension<<=16;
  gap_open_and_extension<<=16;
  dpstate_t *dp[3];
  dp[0]=dp_buffer;
  dp[1]=dp[0]+m+5;
  dp[2]=dp[1]+m+5;
  dp[0][0].idx=dp[1][0].idx=-1;
  int size0=0;
  for (int i=0;i<=n;i++)
    {
      dpstate_t *q0=dp[i&1],*q1=dp[(i+1)&1];
#ifdef DEBUG
      ASSERT(size0>=0 && size0<=MAX_LENGTH);
#endif
      int baseline=est_opt-(n-i);//*(match_score>>16);
      if (baseline<1) baseline=1;
      if (i+7<n && 1>=baseline)
        {
          int pos=1,key=ca[i];
          dpstate_t *q2=dp[2];
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

void Solution5::solve()
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
