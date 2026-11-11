/* Exact integer audit of the CM uniform bound, default mode only.
 * Local character counts are evaluated by their integer formulas.
 * The class-number bound for D>10^7 is rounded upwards using
 * pi>3 and log(2)<7/10. Signed arithmetic overflow traps (-ftrapv).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#define HMAX 10000000
typedef long long ll;
typedef __int128 num;
static int32_t *hcl;
static int nprimes, *plist;
static num ipow(num p, int e) { num r=1; assert(e>=0); while(e--) r*=p; return r; }
static num absnum(num a) { return a<0 ? -a : a; }
static num phi_pk(ll p,int e) { return e==0 ? 1 : ipow(p,e-1)*(p-1); }
static void locinfo(ll p,int tau,int x,int y,num *r,int *v,int *ne) {
  if(tau==0) { *v=x>y?x:y; *ne=x+y; *r=phi_pk(p,x<y?x:y); }
  else if(tau==1) { *v=x; *ne=2*x; *r=x==0?1:ipow(p,x-1)*(p+1); }
  else { *v=(x+1)/2; *ne=x; *r=ipow(p,x/2); }
}
static num Eloc(ll p, int tau, int x, int y, int s, int only2) {
  /* only2: restrict m to divisors of (2) (for the unit correction at p = 2), i.e. rational part v <= 1 and nexp <= 2 */
  num tot = 0;
  for (int dx = 0; dx <= 1; dx++) for (int dy = 0; dy <= (tau == 0 ? 1 : 0); dy++) {
    int mx = x - dx, my = y - dy; if (mx < 0 || my < 0) continue;
    num r; int v, ne; locinfo(p, tau, mx, my, &r, &v, &ne);
    if (only2) { if (!(ne <= 1 || (tau == 0 && mx == 1 && my == 1) || (tau == 1 && mx == 1) || (tau == 2 && mx == 2))) continue; }
    if (s > v) continue;
    int sg = ((dx + dy) % 2) ? -1 : 1;
    tot += sg * r;
  }
  return tot;
}
static num Ehat(ll p, int e, int delta, int tau, int s, int corr) {
  /* corr = 0: sum_f E_p(f) w(f);  corr = 1: the unit-correction local factor */
  num tot = 0; int ymax = (tau == 0) ? e : 0;
  for (int x = 0; x <= e-delta; x++) for (int y = 0; y <= ymax; y++) {
    int nex = tau==0 ? x+y : tau==1 ? 2*x : x;
    if (delta+nex>e) continue;
    num r; int v, ne; locinfo(p, tau, x, y, &r, &v, &ne);
    if (delta + ne > e) continue;
    num w = e - delta - ne + 1;
    num Ev;
    if (!corr) Ev = Eloc(p, tau, x, y, s, 0);
    else if (p == 2) Ev = Eloc(p, tau, x, y, s, 1);
    else { /* odd p: m_p = 1 only: mu(f_p) [s <= 0] */
      int sq = (tau == 0) ? (x <= 1 && y <= 1) : (x <= 1);
      Ev = (sq && s <= 0) ? (((x + y) % 2) ? -1 : 1) : 0; }
    tot += Ev * w;
  }
  return tot;
}

/* local types. odd p: code = a*2 + quad (quad only for a==1). p=2: code = u + 2k. */
static int ntypes(ll p, int e, int *codes) {
  int n = 0;
  if (p == 2) { codes[n++] = 0; if (e >= 2) codes[n++] = 1; for (int k = 1; k <= e - 2; k++) { codes[n++] = 2*k; codes[n++] = 2*k + 1; } }
  else { codes[n++] = 0; codes[n++] = 3; if (p > 3) codes[n++] = 2; for (int a = 2; a <= e; a++) codes[n++] = 2*a; }
  return n;
}
static int s_odd(int delta, int code) { int a = code / 2, quad = code & 1;
  if (delta == 0) return a;
  if (a >= 2) return a;
  if (a == 1) return quad ? 0 : 1;
  return 1; }
static int cond2(int u, int k) { return k == 0 ? 2*u : k + 2; }
static int s_two(int delta, int code, int up) { int u = code & 1, k = code >> 1;
  if (delta == 0) return cond2(u, k);
  if (delta == 2) return k == 0 ? 2*((u+1)%2) : k + 2;
  if (k == 0) return 3;
  if (k == 1) return 2*((u+up)%2);
  return k + 2; }
static ll powmod(ll b, ll e, ll m) { ll r = 1 % m; b %= m; while (e) { if (e&1) r = (__int128)r*b % m; b = (__int128)b*b % m; e >>= 1; } return r; }

static int parmask(ll p, int code) {
  if (p == 2) return (code & 1) ? 2 : 1;
  int a = code / 2, quad = code & 1;
  if (a == 0) return 1;
  if (a == 1 && quad) return (p % 4 == 1) ? 1 : 2;
  if (a == 1) return (p == 5) ? 2 : 3;
  return 3;
}
static int combine(int m1, int m2) { int r = 0;
  for (int x = 0; x < 2; x++) if (m1 >> x & 1) for (int y = 0; y < 2; y++) if (m2 >> y & 1) r |= 1 << (x ^ y);
  return r; }

#define MAXW 16
#define MAXD 8192
#define MAXT 128
static ll P[MAXW]; static int E[MAXW],W,nD;
static num hd[MAXD],L[MAXD][MAXW][MAXT],C[MAXD][MAXW][MAXT];
static num LM[MAXD][MAXW],CM[MAXD][MAXW];
static int TN[MAXW],TC[MAXW][MAXT],TP[MAXW][MAXT],sel[MAXW];
static num threshold; static int found;
static num part[MAXD],partC[MAXD];
static void bb(int k,int parity) {
  if(found) return;
  if(k==W) {
    if(!(parity&1)) return;
    num value=0;
    for(int d=0;d<nD;d++) {
      num a=1,b=1;
      for(int j=0;j<W;j++) { a*=L[d][j][sel[j]]; b*=C[d][j][sel[j]]; }
      value+=hd[d]*(a-b);
    }
    if(2*value>=threshold) found=1;
    return;
  }
  for(int t=0;t<TN[k];t++) {
    num save[MAXD],saveC[MAXD],opt=0;
    for(int d=0;d<nD;d++) {
      save[d]=part[d]; saveC[d]=partC[d];
      part[d]*=L[d][k][t]; partC[d]*=absnum(C[d][k][t]);
      num a=part[d],b=partC[d];
      for(int j=k+1;j<W;j++) { a*=LM[d][j]; b*=CM[d][j]; }
      opt+=hd[d]*(a+b);
    }
    sel[k]=t;
    if(2*opt>=threshold) bb(k+1,combine(parity,TP[k][t]));
    for(int d=0;d<nD;d++) { part[d]=save[d]; partC[d]=saveC[d]; }
    if(found) return;
  }
}
static num classbound(ll D) {
  if(D<=HMAX) return hcl[D];
  ll lo=0,hi=1;
  while((num)hi*hi<D) hi*=2;
  while(hi-lo>1) { ll m=(hi+lo)/2; if((num)m*m>=D) hi=m; else lo=m; }
  int bits=0; for(ll d=D;d;d/=2) bits++;
  return ((num)hi*(7*bits+20))/30+1;
}
int main(void) {
  hcl=calloc(HMAX+1,sizeof(int32_t)); assert(hcl);
  for(ll a=1;3*a*a<=HMAX;a++) for(ll b=-a+1;b<=a;b++)
    for(ll c=a;4*a*c-b*b<=HMAX;c++) {
      if(c==a && b<0) continue;
      ll d=4*a*c-b*b; if(d>=3) hcl[d]++;
    }
  fprintf(stderr,"Exact class numbers initialized.\n");
  plist=malloc(100000*sizeof(int)); nprimes=0;
  char *sv=calloc(1000001,1); assert(plist && sv);
  for(int i=2;i<=1000000;i++) if(!sv[i]) {
    plist[nprimes++]=i;
    for(ll j=(ll)i*i;j<=1000000;j+=i) sv[j]=1;
  }
  free(sv);
  ll N; long checked=0,fail=0;
  while(scanf("%lld",&N)==1) {
    assert(N>=1 && N<=266378112000LL); checked++;
    ll n=N; W=0;
    for(int i=0;i<nprimes && (ll)plist[i]*plist[i]<=n;i++) {
      ll p=plist[i]; if(n%p) continue;
      assert(W<MAXW); P[W]=p; E[W]=0;
      while(n%p==0) { n/=p; E[W]++; } W++;
    }
    if(n>1) { assert(W<MAXW); P[W]=n; E[W++]=1; }
    num psi=N,nu2=1,nu3=1,nuinf=1;
    int four=0,nine=0;
    for(int k=0;k<W;k++) {
      ll p=P[k]; int e=E[k]; psi=psi/p*(p+1);
      if(p==2 && e>=2) four=1;
      if(p==3 && e>=2) nine=1;
      if(p!=2) nu2*=p%4==1?2:0;
      if(p!=3) nu3*=p%3==1?2:0;
      num ss=0;
      for(int j=0;j<=e;j++) ss+=phi_pk(p,j<e-j?j:e-j);
      nuinf*=ss;
      TN[k]=ntypes(p,e,TC[k]); assert(TN[k]<=MAXT);
      for(int t=0;t<TN[k];t++) TP[k][t]=parmask(p,TC[k][t]);
    }
    if(four) nu2=0;
    if(nine) nu3=0;
    num numerator=psi-3*nu2-4*nu3-6*nuinf;
    assert(numerator%12==0); threshold=numerator/12;
    nD=0; int oi[MAXW],no=0,e2=0;
    for(int k=0;k<W;k++) if(P[k]==2) e2=E[k]; else oi[no++]=k;
    for(int mask=0;mask<(1<<no);mask++) {
      ll d=1; int n3=0;
      for(int t=0;t<no;t++) if(mask>>t&1) { d*=P[oi[t]]; if(P[oi[t]]%4==3) n3++; }
      ll cand[3]; int nc=0;
      if(d%4==3) cand[nc++]=d;
      if(e2>=2 && d%4==1) cand[nc++]=4*d;
      if(e2>=3) cand[nc++]=8*d;
      for(int cc=0;cc<nc;cc++) {
        ll D=cand[cc]; assert(nD<MAXD); hd[nD]=classbound(D);
        int up=n3%2==0?1:0;
        for(int k=0;k<W;k++) {
          ll p=P[k]; int e=E[k],delta,tau;
          if(p==2) { delta=D%8==0?3:D%4==0?2:0; tau=delta?2:((8-D%8)%8==1?0:1); }
          else { delta=D%p==0; tau=delta?2:(powmod((p-D%p)%p,(p-1)/2,p)==1?0:1); }
          LM[nD][k]=0; CM[nD][k]=0;
          for(int t=0;t<TN[k];t++) {
            int code=TC[k][t]; int ss=p==2?s_two(delta,code,up):s_odd(delta,code);
            num v=Ehat(p,e,delta,tau,ss,0);
            num c=D>4?Ehat(p,e,delta,tau,ss,1):0;
            assert(v>=0);
            L[nD][k][t]=v; C[nD][k][t]=c;
            if(v>LM[nD][k]) LM[nD][k]=v;
            if(absnum(c)>CM[nD][k]) CM[nD][k]=absnum(c);
          }
        }
        nD++;
      }
    }
    num opt=0;
    for(int d=0;d<nD;d++) {
      num a=hd[d],b=hd[d];
      for(int k=0;k<W;k++) { a*=LM[d][k]; b*=CM[d][k]; }
      opt+=a+b; part[d]=1; partC[d]=1;
    }
    found=threshold<=0;
    if(!found && 2*opt>=threshold) bb(0,1);
    if(found) { fail++; printf("FAIL N=%lld\n",N); fflush(stdout); }
    if(checked%100000==0) { fprintf(stderr,"checked=%ld last_N=%lld\n",checked,N); }
  }
  printf("SUMMARY checked=%ld fail=%ld\n",checked,fail);
  return 0;
}
