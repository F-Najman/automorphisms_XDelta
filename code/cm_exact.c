/* Exact check of (P_N): 2 M(N) < g_0(N) - 1, general N (Section 6 of the paper), for N read from stdin.
 * M(N) <= max over local types t = (t_p) of chi of
 *   V(t) = sum_{D | N fundamental} h(-D) [ prod_p L_p(D, t_p) - Excl_D(t) ],
 * L_p(D,t) = sum over (i,j) with delta + 2i + j <= e, i + j >= s_p(D,t), of rho(p^i) a(p^j) (e - delta - 2i - j + 1),
 * with the true splitting type of p in Q(sqrt -D); Excl_D(t) = the terms with conductor dividing (2) (no Hecke
 * character of infinity type alpha -> alpha exists there, since -1 = 1 mod f).
 * h(-D): exact (reduced forms) for D <= HMAX, else floor(sqrt(D)(log D + 2)/pi) + 1.
 * Branch and bound over types with the optimistic bound sum_D h prod(assigned) prod(max over unassigned) (no Excl).
 * Output per failing N: "FAIL N=... V=... theta=...", and a summary.   Build: gcc -O2 -fopenmp -o cm_exact cm_exact.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#define HMAX 10000000
static int32_t *hcl;
typedef long long ll;
static int nprimes; static int *plist;

static double rho(ll p, int i, int tau) { /* tau: 0 split, 1 inert, 2 ram */
  if (i == 0) return 1.0; double q = pow((double)p, i);
  return tau == 0 ? q * (1 - 1.0/p) : tau == 1 ? q * (1 + 1.0/p) : q;
}
static int amult(int j, int tau) { if (j == 0) return 1; if (tau == 0) return 2; if (tau == 2) return j == 1; return 0; }
static double Lf(ll p, int e, int delta, int tau, int s, int excl_only) {
  /* excl_only: only the (i,j) with conductor part dividing (2) (used at p = 2) */
  double tot = 0;
  for (int i = 0; 2*i <= e; i++) for (int j = 0; delta + 2*i + j <= e; j++) {
    if (i + j < s) continue;
    int am = amult(j, tau); if (!am) continue;
    if (excl_only && !(i + j <= 1)) continue;          /* 1, primes above 2, and (2) = c=2 */
    tot += rho(p, i, tau) * am * (e - delta - 2*i - j + 1);
  }
  return tot;
}

/* ---- exact-conductor (Moebius) local factors, valid for D > 4 (w = 2) ----
 * local ideals f_p: split p^x pbar^y, inert (p)^x, ram q^x.  norm exponent nexp, rational part v = v_p(f_p cap Z).
 * r(m) = phi(m)/phi(p^v(m)); S_loc(m) nonempty iff s <= v(m).  E_p(f) = sum_{m | f, f/m squarefree} mu(f/m) r(m) [s <= v(m)]. */
static double phi_pk(double p, int k) { return k == 0 ? 1.0 : pow(p, k) - pow(p, k - 1); }
static void locinfo(ll p, int tau, int x, int y, double *r, int *v, int *nexp) {
  if (tau == 0) { *v = x > y ? x : y; *nexp = x + y; *r = phi_pk(p, x) * phi_pk(p, y) / phi_pk(p, *v); }
  else if (tau == 1) { *v = x; *nexp = 2 * x; *r = (x == 0) ? 1.0 : (pow(p, 2*x) - pow(p, 2*x - 2)) / phi_pk(p, x); }
  else { *v = (x + 1) / 2; *nexp = x; *r = phi_pk(p, x) / phi_pk(p, *v); }
}
static double Eloc(ll p, int tau, int x, int y, int s, int only2) {
  /* only2: restrict m to divisors of (2) (for the unit correction at p = 2), i.e. rational part v <= 1 and nexp <= 2 */
  double tot = 0;
  for (int dx = 0; dx <= 1; dx++) for (int dy = 0; dy <= (tau == 0 ? 1 : 0); dy++) {
    int mx = x - dx, my = y - dy; if (mx < 0 || my < 0) continue;
    double r; int v, ne; locinfo(p, tau, mx, my, &r, &v, &ne);
    if (only2) { if (!(ne <= 1 || (tau == 0 && mx == 1 && my == 1) || (tau == 1 && mx == 1) || (tau == 2 && mx == 2))) continue; }
    if (s > v) continue;
    int sg = ((dx + dy) % 2) ? -1 : 1;
    tot += sg * r;
  }
  return tot;
}
static double Ehat(ll p, int e, int delta, int tau, int s, int corr) {
  /* corr = 0: sum_f E_p(f) w(f);  corr = 1: the unit-correction local factor */
  double tot = 0; int ymax = (tau == 0) ? e : 0;
  for (int x = 0; x <= 2 * e; x++) for (int y = 0; y <= ymax; y++) {
    double r; int v, ne; locinfo(p, tau, x, y, &r, &v, &ne);
    if (delta + ne > e) continue;
    double w = e - delta - ne + 1;
    double Ev;
    if (!corr) Ev = Eloc(p, tau, x, y, s, 0);
    else if (p == 2) Ev = Eloc(p, tau, x, y, s, 1);
    else { /* odd p: m_p = 1 only: mu(f_p) [s <= 0] */
      int sq = (tau == 0) ? (x <= 1 && y <= 1) : (x <= 1);
      Ev = (sq && s <= 0) ? (((x + y) % 2) ? -1.0 : 1.0) : 0.0; }
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
  if (delta == 0) return a; if (a >= 2) return a; if (a == 1) return quad ? 0 : 1; return 1; }
static int cond2(int u, int k) { return k == 0 ? 2*u : k + 2; }
static int s_two(int delta, int code, int up) { int u = code & 1, k = code >> 1;
  if (delta == 0) return cond2(u, k);
  if (delta == 2) return k == 0 ? 2*((u+1)%2) : k + 2;
  if (k == 0) return 3; if (k == 1) return 2*((u+up)%2); return k + 2; }
static ll powmod(ll b, ll e, ll m) { ll r = 1 % m; b %= m; while (e) { if (e&1) r = (__int128)r*b % m; b = (__int128)b*b % m; e >>= 1; } return r; }

#define MAXW 16
#define MAXD 8192
#define MAXT 128
static ll P[MAXW]; static int E[MAXW]; static int W;
static int nD; static ll Dv[MAXD]; static double hD[MAXD];
static double Lt[MAXD][MAXW][MAXT];   /* L_p(D, type) */
static double Lmax[MAXD][MAXW]; static double Cmax[MAXD][MAXW];
static int Tn[MAXW]; static int Tc[MAXW][MAXT]; static int Tpar[MAXW][MAXT]; static double Tlam[MAXW][MAXT]; static double lammax[MAXW];
static int dimmode = 0; static double psi12, c34;
static double ex_odd[MAXD][MAXW][MAXT]; /* odd p: (e-delta+1) if s==0 else 0 ; p=2: excl sum */
static double theta; static int found; static double bestV; static int maxmode = 0;
static int tsel[MAXW];
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
static int parstack[MAXW + 1];
static double lamstack[MAXW + 1];
static double lamCO(ll p, int r, int s2) { /* Cohen-Oesterle lambda(r, s, p) */
  if (2*s2 <= r) { if (r % 2 == 0) { int rr = r/2; return pow((double)p, rr) + (rr ? pow((double)p, rr-1) : 0); } else return 2*pow((double)p, r/2); }
  return 2*pow((double)p, r - s2); }

static double evalV(void) {
  double V = 0;
  for (int d = 0; d < nD; d++) {
    double pr = 1, ex = 1;
    for (int k = 0; k < W; k++) { pr *= Lt[d][k][tsel[k]]; ex *= ex_odd[d][k][tsel[k]]; }
    V += hD[d] * (pr - ex);
  }
  return V;
}
static double part[MAXD]; static double partC[MAXD];
static void bb(int k) {
  if (found) return;
  if (k == W) { if (!(parstack[W] & 1)) return; double V = evalV();
    if (dimmode) { double dl = psi12 - 0.5 * lamstack[W] - c34; double val = 2*V - (dl - 1); if (val > bestV) bestV = val; if (val > 0) found = 1; return; }
    if (V > bestV) bestV = V; if (!maxmode && V >= theta) found = 1; return; }
  for (int ti = 0; ti < Tn[k]; ti++) {
    double save[MAXD], saveC[MAXD]; double opt = 0;
    for (int d = 0; d < nD; d++) { save[d] = part[d]; saveC[d] = partC[d];
      part[d] *= Lt[d][k][ti]; partC[d] *= fabs(ex_odd[d][k][ti]);
      double o = part[d], oc = partC[d]; for (int m = k+1; m < W; m++) { o *= Lmax[d][m]; oc *= Cmax[d][m]; } opt += hD[d] * (o + oc); }
    tsel[k] = ti; parstack[k+1] = combine(parstack[k], Tpar[k][ti]); lamstack[k+1] = lamstack[k] * Tlam[k][ti];
    int go;
    if (dimmode) { double lr = lamstack[k+1]; for (int m = k+1; m < W; m++) lr *= lammax[m];
      go = (2*opt - (psi12 - 0.5*lr - c34 - 1) > 0); }
    else go = maxmode ? (opt > bestV) : (opt >= theta);
    if (go) bb(k + 1);
    for (int d = 0; d < nD; d++) { part[d] = save[d]; partC[d] = saveC[d]; }
    if (found) return;
  }
}

int main(int argc, char **argv) {
  if (argc > 1 && !strcmp(argv[1], "max")) maxmode = 1;
  if (argc > 1 && !strcmp(argv[1], "dim")) dimmode = 1;
  hcl = calloc(HMAX + 1, sizeof(int32_t));
  for (ll a = 1; 3*a*a <= HMAX; a++) for (ll b = -a + 1; b <= a; b++) for (ll c = a; 4*a*c - b*b <= HMAX; c++) {
    if (c == a && b < 0) continue; ll d = 4*a*c - b*b; if (d >= 3) hcl[d]++; }
  fprintf(stderr, "class numbers done: h(23)=%d h(999983)=%d\n", hcl[23], hcl[999983]);
  plist = malloc(sizeof(int) * 600000); nprimes = 0;
  { char *sv = calloc(1000001, 1); for (int i = 2; i <= 1000000; i++) if (!sv[i]) { plist[nprimes++] = i; for (ll j = (ll)i*i; j <= 1000000; j += i) sv[j] = 1; } free(sv); }
  ll N; long ncheck = 0, nfail = 0;
  while (scanf("%lld", &N) == 1) {
    ncheck++;
    /* factor */
    ll n = N; W = 0;
    for (int i = 0; i < nprimes && (ll)plist[i]*plist[i] <= n; i++) { ll p = plist[i]; if (n % p == 0) { P[W] = p; E[W] = 0; while (n % p == 0) { n /= p; E[W]++; } W++; } }
    if (n > 1) { P[W] = n; E[W] = 1; W++; }
    /* sort primes so that 2 comes first (already ascending) */
    /* g_0 - 1 */
    double psi = N; ll nu2 = 1, nu3 = 1, nuinf = 1; int four = 0, nine = 0;
    for (int k = 0; k < W; k++) { ll p = P[k]; int e = E[k]; psi *= (1.0 + 1.0/p);
      if (p == 2 && e >= 2) four = 1; if (p == 3 && e >= 2) nine = 1;
      if (p != 2) { ll kr = (p % 4 == 1) ? 1 : -1; nu2 *= 1 + kr; }
      if (p != 3) { ll kr = (p % 3 == 1) ? 1 : (p == 3 ? 0 : -1); nu3 *= 1 + kr; }
      ll s = 0; for (int j = 0; j <= e; j++) { int m = j < e - j ? j : e - j; s += m == 0 ? 1 : (ll)(pow((double)p, m) - pow((double)p, m - 1)); } nuinf *= s; }
    if (four) nu2 = 0; if (nine) nu3 = 0;
    double g0m1 = psi/12.0 - nu2/4.0 - nu3/3.0 - nuinf/2.0;
    theta = g0m1 / 2.0;
    psi12 = psi / 12.0;
    { ll s4 = four ? 0 : 1, s3 = nine ? 0 : 1;   /* #sol of x^2+1 and x^2+x+1 mod N */
      for (int k = 0; k < W; k++) { ll p = P[k];
        if (p == 2) { s4 *= (E[k] == 1) ? 1 : 0; s3 *= 0; }
        else if (p == 3) { s4 *= 0; s3 *= (E[k] == 1) ? 1 : 0; }
        else { s4 *= (p % 4 == 1) ? 2 : 0; s3 *= (p % 3 == 1) ? 2 : 0; } }
      c34 = s4 / 4.0 + s3 / 3.0; }
    /* fundamental discriminants D | N */
    nD = 0; int oddidx[MAXW], no = 0, e2 = 0;
    for (int k = 0; k < W; k++) { if (P[k] == 2) e2 = E[k]; else oddidx[no++] = k; }
    for (int mask = 0; mask < (1 << no); mask++) {
      ll d = 1; int n3 = 0;
      for (int t = 0; t < no; t++) if (mask >> t & 1) { d *= P[oddidx[t]]; if (P[oddidx[t]] % 4 == 3) n3++; }
      ll cand[3]; int nc = 0;
      if (d % 4 == 3) cand[nc++] = d;
      if (e2 >= 2 && d % 4 == 1) cand[nc++] = 4*d;
      if (e2 >= 3) cand[nc++] = 8*d;
      for (int c = 0; c < nc; c++) { ll D = cand[c]; if (D < 3) continue;
        if (nD >= MAXD) { fprintf(stderr, "MAXD\n"); exit(1); }
        Dv[nD] = D;
        hD[nD] = D <= HMAX ? hcl[D] : floor(sqrt((double)D) * (log((double)D) + 2) / M_PI) + 1;
        int up = (n3 % 2 == 0) ? 1 : 0;   /* eps_{F,2}(-1) = -(-1)^{n3}: +1 -> chi_8 (up=0), -1 -> chi_-8 (up=1) */
        for (int k = 0; k < W; k++) {
          ll p = P[k]; int e = E[k]; int delta, tau;
          if (p == 2) { delta = (D % 8 == 0) ? 3 : (D % 4 == 0) ? 2 : 0; if (delta) tau = 2; else tau = ((8 - D % 8) % 8 == 1) ? 0 : 1; }
          else { delta = (D % p == 0); if (delta) tau = 2; else { ll x = (p - D % p) % p; tau = powmod(x, (p-1)/2, p) == 1 ? 0 : 1; } }
          Tn[k] = ntypes(p, e, Tc[k]);
          lammax[k] = 0;
          for (int ti = 0; ti < Tn[k]; ti++) { Tpar[k][ti] = parmask(p, Tc[k][ti]);
            int code = Tc[k][ti]; int a = (p == 2) ? cond2(code & 1, code >> 1) : code / 2;
            Tlam[k][ti] = lamCO(p, e, a); if (Tlam[k][ti] > lammax[k]) lammax[k] = Tlam[k][ti]; }
          if (Tn[k] > MAXT) { fprintf(stderr, "MAXT\n"); exit(1); }
          double mx = 0, cmx = 0;
          for (int ti = 0; ti < Tn[k]; ti++) {
            int code = Tc[k][ti];
            int s = (p == 2) ? s_two(delta, code, up) : s_odd(delta, code);
            double v;
            if (D > 4) { v = Ehat(p, e, delta, tau, s, 0); ex_odd[nD][k][ti] = Ehat(p, e, delta, tau, s, 1); }
            else { /* D = 3, 4 (h = 1): exact-conductor count of characters of (O/f)^x with prescribed
                      restriction to (Z/f)^x, dropping the extra unit conditions (an upper bound) */
              v = Ehat(p, e, delta, tau, s, 0); ex_odd[nD][k][ti] = 0.0; }
            Lt[nD][k][ti] = v; if (v > mx) mx = v; if (fabs(ex_odd[nD][k][ti]) > cmx) cmx = fabs(ex_odd[nD][k][ti]);
          }
          Lmax[nD][k] = mx; Cmax[nD][k] = cmx;
        }
        nD++;
      }
    }
    found = 0; bestV = 0;
    double opt = 0; for (int d = 0; d < nD; d++) { double o = hD[d], oc = hD[d]; for (int k = 0; k < W; k++) { o *= Lmax[d][k]; oc *= Cmax[d][k]; } opt += o + oc; }
    parstack[0] = 1; lamstack[0] = 1.0;
    if (maxmode) { for (int d = 0; d < nD; d++) { part[d] = 1; partC[d] = 1; } bb(0); printf("MAX N=%lld Mplus=%.2f theta=%.2f\n", N, bestV, theta); continue; }
    if (dimmode) { for (int d = 0; d < nD; d++) { part[d] = 1; partC[d] = 1; } bb(0); }
    else if (theta <= 0) found = 1;           /* g_0 <= 1: (P_N) cannot hold */
    else if (opt >= theta) { for (int d = 0; d < nD; d++) { part[d] = 1; partC[d] = 1; } bb(0); }
    if (found) { nfail++; printf("FAIL N=%lld theta=%.2f bestV>=%.2f opt=%.2f nD=%d\n", N, theta, bestV, opt, nD); fflush(stdout); }
  }
  printf("SUMMARY checked=%ld fail=%ld\n", ncheck, nfail);
  return 0;
}
