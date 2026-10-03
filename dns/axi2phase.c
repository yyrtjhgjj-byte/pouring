/*
 * axi2phase.c — 軸対称・非圧縮・二相流（水／空気）の直接数値シミュレーション
 *
 * ウォーターサーバーの水柱がコップの水面に突っ込む瞬間を、表面張力込みで解く。
 *
 *   支配方程式   ρ(∂u/∂t + u·∇u) = -∇p + ∇·(2μD) + σκδ_s n + ρg,   ∇·u = 0
 *   座標         (r, z) 軸対称（r=0 が対称軸）、z 上向き、重力 -z
 *   格子         一様 MAC スタガード格子（u: r 面, w: z 面, p・c: セル中心）
 *   界面         PLIC-VOF（Youngs 法線）+ Weymouth & Yue (2010) の質量保存型方向分割移流
 *   曲率         高さ関数法（軸対称項込み）、だめなセルは近傍平均 → 平滑化 CSF にフォールバック
 *   表面張力     balanced-force CSF（圧力勾配と同じ面で離散化 → 静止液滴で寄生流がほぼ出ない）
 *   移流         MC 制限付き 2 次風上 + SSP-RK2
 *   粘性         陽的、変粘性の応力発散形（軸対称の -2μu/r² 項込み）
 *   圧力         変係数ポアソン（密度比 ~830）を、Galerkin 幾何マルチグリッド前処理付き
 *                Flexible CG で解く。OpenMP 並列。
 *
 * 境界: r=0 対称軸 / r=Lr 滑り壁（コップの側壁）/ z=0 滑り壁（底）/
 *       z=Lz 上端は開放（p=0）で、r<a_in から水柱が流入（流入停止時刻を指定可）
 *
 * ビルド: gcc -O3 -march=native -fopenmp -o axi2phase axi2phase.c -lm
 * 実行:   ./axi2phase params.txt
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define NG 4
#define SQ(x) ((x) * (x))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define CLAMP(x, a, b) MIN(MAX(x, a), b)

/* ------------------------------------------------------------------ パラメータ */
typedef struct {
  int nr, nz;
  double Lr;
  /* 物性: 1 = 水, 2 = 空気 */
  double rho1, rho2, mu1, mu2, sigma, grav;
  /* 初期条件 */
  double pool_depth;           /* 水面の高さ（<=0 でプールなし） */
  double jet_radius;           /* 流入水柱の半径 */
  double jet_velocity;         /* 流入速度（下向き正） */
  double jet_head_radius;      /* 先頭の水塊の半径（0 なら水柱半径） */
  double jet_tip_z;            /* 初期の水柱先端の高さ（これより上が水柱） */
  double jet_stop_time;        /* 流入を止める時刻（負なら止めない） */
  int jet_init;                /* 1: 初期に水柱をドメイン内に置く */
  int jet_profile;             /* 1: 初期水柱を重力加速の Bernoulli 形状にする（細りながら加速）*/
  double drop_radius, drop_z, drop_velocity;   /* 単独液滴（検証用） */
  double drop_p2;                               /* 液滴の P2 変形（振動液滴テスト）*/
  double bubble_radius, bubble_z;              /* 単独気泡（検証用） */
  double cavity_radius, cavity_depth;          /* 水面のくぼみ（気泡破裂の初期形状） */
  /* 数値 */
  double t_end, cfl, out_dt, tol;
  int max_steps, mg_pre, mg_post, verbose;
  char out_dir[512];
} Params;

static void params_default(Params *P) {
  memset(P, 0, sizeof(*P));
  P->nr = 128; P->nz = 256; P->Lr = 0.032;
  P->rho1 = 998.2; P->rho2 = 1.204; P->mu1 = 1.002e-3; P->mu2 = 1.81e-5;
  P->sigma = 0.0727; P->grav = 9.80665;
  P->pool_depth = 0.03;
  P->jet_radius = 0.0; P->jet_velocity = 0.0; P->jet_head_radius = 0.0;
  P->jet_tip_z = 0.0; P->jet_stop_time = -1.0; P->jet_init = 1;
  P->t_end = 0.1; P->cfl = 0.4; P->out_dt = 1e-3; P->tol = 1e-5;
  P->max_steps = 10000000; P->mg_pre = 2; P->mg_post = 2; P->verbose = 1;
  strcpy(P->out_dir, "out");
}

static void params_read(Params *P, const char *fn) {
  FILE *f = fopen(fn, "r");
  if (!f) { perror(fn); exit(1); }
  char key[128], val[512], line[1024];
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '#' || sscanf(line, "%127s %511s", key, val) != 2) continue;
#define PD(name) else if (!strcmp(key, #name)) P->name = atof(val);
#define PI(name) else if (!strcmp(key, #name)) P->name = atoi(val);
    if (0) {}
    PI(nr) PI(nz) PD(Lr) PD(rho1) PD(rho2) PD(mu1) PD(mu2) PD(sigma) PD(grav)
    PD(pool_depth) PD(jet_radius) PD(jet_velocity) PD(jet_head_radius) PD(jet_tip_z)
    PD(jet_stop_time) PI(jet_init) PI(jet_profile) PD(drop_radius) PD(drop_z) PD(drop_velocity) PD(drop_p2)
    PD(bubble_radius) PD(bubble_z) PD(cavity_radius) PD(cavity_depth)
    PD(t_end) PD(cfl) PD(out_dt) PD(tol) PI(max_steps) PI(mg_pre) PI(mg_post) PI(verbose)
    else if (!strcmp(key, "out_dir")) strncpy(P->out_dir, val, sizeof(P->out_dir) - 1);
    else fprintf(stderr, "warning: unknown key %s\n", key);
  }
  fclose(f);
}

/* ------------------------------------------------------------------ 格子と配列 */
static int NR, NZ, NI, NJ;
static double DX;
#define ID(i, j) (((i) + NG) * NJ + (j) + NG)

static double *alloc_field(void) {
  double *a = calloc((size_t)NI * NJ, sizeof(double));
  if (!a) { fprintf(stderr, "out of memory\n"); exit(1); }
  return a;
}
static inline double rc(int i) { return (i + 0.5) * DX; } /* セル中心の r */
static inline double rf(int i) { return i * DX; }         /* i 番目の r 面（セル i の左）*/

static Params P;
static double t_now = 0.0;
static double T_vof = 0, T_curv = 0, T_mom = 0, T_proj = 0;
static double wclock(void) {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / CLOCKS_PER_SEC;
#endif
}
/* 主要な場 */
static double *c, *u, *w, *p, *kap, *mx_, *my_, *alph;
static double *u_s, *w_s, *u1, *w1, *du, *dw;          /* 作業配列 */
static double *rho_c, *mu_c, *csm, *csm2;
static int *hf_ok;

/* 流入している r 範囲か */
static int jet_on(void) {
  return P.jet_radius > 0 && P.jet_velocity > 0 && (P.jet_stop_time < 0 || t_now < P.jet_stop_time);
}
/* 流入面 i が水柱にかかる割合（軸対称の体積重み） */
static double inflow_frac(int i) {
  double a = P.jet_radius, r0 = rf(i), r1 = rf(i + 1);
  if (a <= r0) return 0.0;
  if (a >= r1) return 1.0;
  return (a * a - r0 * r0) / (r1 * r1 - r0 * r0);
}

/* ------------------------------------------------------------------ 境界条件 */
static void bc_c(double *f) {
  /* 軸と側壁と底は鏡像、上端は流入部 = 1、他は勾配ゼロ */
#pragma omp parallel for
  for (int j = -NG; j < NZ + NG; j++)
    for (int g = 1; g <= NG; g++) {
      f[ID(-g, j)] = f[ID(g - 1, j)];
      f[ID(NR - 1 + g, j)] = f[ID(NR - g, j)];
    }
  int on = jet_on();
#pragma omp parallel for
  for (int i = -NG; i < NR + NG; i++) {
    int ii = CLAMP(i, 0, NR - 1);
    for (int g = 1; g <= NG; g++) {
      f[ID(i, -g)] = f[ID(i, g - 1)];
      double top = f[ID(i, NZ - 1)];
      if (on) {
        double fr = inflow_frac(ii);
        if (fr > 0) top = fr;
      }
      f[ID(i, NZ - 1 + g)] = top;
    }
  }
}

static void bc_vel(double *U, double *W) {
  int on = jet_on();
  /* r 方向: 軸で u 奇, w 偶; 側壁（滑り）で u 奇, w 偶 */
#pragma omp parallel for
  for (int j = -NG; j < NZ + NG; j++) {
    U[ID(0, j)] = 0.0;
    U[ID(NR, j)] = 0.0;
    for (int g = 1; g <= NG; g++) {
      U[ID(-g, j)] = -U[ID(g, j)];
      U[ID(NR + g, j)] = -U[ID(NR - g, j)];
      W[ID(-g, j)] = W[ID(g - 1, j)];
      W[ID(NR - 1 + g, j)] = W[ID(NR - g, j)];
    }
  }
  /* z 方向: 底で w 奇, u 偶; 上端は流入（w 固定, u=0）か流出（勾配ゼロ）*/
#pragma omp parallel for
  for (int i = -NG; i < NR + NG + 1; i++) {
    int ii = CLAMP(i, 0, NR - 1);
    W[ID(i, 0)] = 0.0;
    for (int g = 1; g <= NG; g++) {
      W[ID(i, -g)] = -W[ID(i, g)];
      U[ID(i, -g)] = U[ID(i, g - 1)];
    }
    double fr = on ? inflow_frac(ii) : 0.0;
    if (fr > 0) {
      W[ID(i, NZ)] = -P.jet_velocity;
      for (int g = 1; g <= NG; g++) {
        W[ID(i, NZ + g)] = W[ID(i, NZ)];
        U[ID(i, NZ - 1 + g)] = -U[ID(i, NZ - g)];
      }
    } else {
      for (int g = 1; g <= NG; g++) {
        W[ID(i, NZ + g)] = W[ID(i, NZ)];
        U[ID(i, NZ - 1 + g)] = U[ID(i, NZ - 1)];
      }
    }
  }
}

/* 上端の面 i が流入面か */
static inline int is_inflow_face(int i) { return jet_on() && inflow_frac(i) > 0; }

/* ------------------------------------------------------------------ PLIC 幾何 */
/* 単位正方形 [-.5,.5]^2 で n·x <= alpha の面積（Basilisk の line_area と同じ規約）*/
static double line_area(double nx, double ny, double alpha) {
  alpha += (nx + ny) / 2.;
  if (nx < 0.) { alpha -= nx; nx = -nx; }
  if (ny < 0.) { alpha -= ny; ny = -ny; }
  if (alpha <= 0.) return 0.;
  if (alpha >= nx + ny) return 1.;
  double area;
  if (nx < 1e-10) area = alpha / ny;
  else if (ny < 1e-10) area = alpha / nx;
  else {
    double v = alpha * alpha, a = alpha - nx;
    if (a > 0.) v -= a * a;
    a = alpha - ny;
    if (a > 0.) v -= a * a;
    area = v / (2. * nx * ny);
  }
  return CLAMP(area, 0., 1.);
}
/* 体積率 cf を与える直線定数 alpha（|nx|+|ny|=1） */
static double line_alpha(double cf, double nx, double ny) {
  double n1 = fabs(nx), n2 = fabs(ny), alpha;
  if (n1 > n2) { double t = n1; n1 = n2; n2 = t; }
  cf = CLAMP(cf, 0., 1.);
  double v1 = n1 / 2.;
  if (cf <= v1 / n2) alpha = sqrt(2. * cf * n1 * n2);
  else if (cf <= 1. - v1 / n2) alpha = cf * n2 + v1;
  else alpha = n1 + n2 - sqrt(2. * n1 * n2 * (1. - cf));
  if (nx < 0.) alpha += nx;
  if (ny < 0.) alpha += ny;
  return alpha - (nx + ny) / 2.;
}
/* セル座標の矩形 [x0,x1]x[y0,y1] の中の液体割合 */
static double rect_fraction(double nx, double ny, double alpha, double x0, double x1, double y0, double y1) {
  alpha -= nx * (x0 + x1) / 2. + ny * (y0 + y1) / 2.;
  return line_area(nx * (x1 - x0), ny * (y1 - y0), alpha);
}

/* Youngs 法線（液体の外向き = -∇c, L1 正規化）と alpha を作る */
static void reconstruct(void) {
  bc_c(c);
#pragma omp parallel for
  for (int i = -1; i < NR + 1; i++)
    for (int j = -1; j < NZ + 1; j++) {
      int k = ID(i, j);
      double cc = c[k];
      if (cc <= 0. || cc >= 1.) { mx_[k] = 0; my_[k] = 0; alph[k] = 0; continue; }
      double gx = (c[ID(i + 1, j + 1)] + 2 * c[ID(i + 1, j)] + c[ID(i + 1, j - 1)])
                - (c[ID(i - 1, j + 1)] + 2 * c[ID(i - 1, j)] + c[ID(i - 1, j - 1)]);
      double gy = (c[ID(i + 1, j + 1)] + 2 * c[ID(i, j + 1)] + c[ID(i - 1, j + 1)])
                - (c[ID(i + 1, j - 1)] + 2 * c[ID(i, j - 1)] + c[ID(i - 1, j - 1)]);
      double nx = -gx, ny = -gy, nn = fabs(nx) + fabs(ny);
      if (nn < 1e-12) { nx = 1.; ny = 0.; nn = 1.; }
      nx /= nn; ny /= nn;
      mx_[k] = nx; my_[k] = ny;
      alph[k] = line_alpha(cc, nx, ny);
    }
}

/* ------------------------------------------------------------------ VOF 移流 */
static double mass_lost = 0.0;

static double vof_flux_frac(int k, double s, int dir) {
  /* donor セル k、掃引割合 s(>0: 正方向へ出る) */
  double cc = c[k];
  if (cc <= 0. || cc >= 1.) return cc;
  double nx = mx_[k], ny = my_[k], a = alph[k], as = fabs(s);
  if (dir == 0) {
    return s > 0 ? rect_fraction(nx, ny, a, 0.5 - as, 0.5, -0.5, 0.5)
                 : rect_fraction(nx, ny, a, -0.5, -0.5 + as, -0.5, 0.5);
  } else {
    return s > 0 ? rect_fraction(nx, ny, a, -0.5, 0.5, 0.5 - as, 0.5)
                 : rect_fraction(nx, ny, a, -0.5, 0.5, -0.5, -0.5 + as);
  }
}

static void vof_sweep(int dir, double dt, const double *cc0) {
  reconstruct();
  /* 面フラックス（u_s/w_s を流用） */
  double *F = du;
  if (dir == 0) {
#pragma omp parallel for
    for (int i = 0; i <= NR; i++)
      for (int j = 0; j < NZ; j++) {
        double uf = u[ID(i, j)], s = uf * dt / DX;
        int donor = s > 0 ? ID(i - 1, j) : ID(i, j);
        F[ID(i, j)] = (uf == 0.) ? 0. : vof_flux_frac(donor, s, 0) * uf;
      }
    double ml = 0.0;
#pragma omp parallel for reduction(+ : ml)
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = ID(i, j);
        double r0 = rf(i), r1 = rf(i + 1), r = rc(i);
        double cn = c[k] + dt / (DX * r) * ((r0 * F[ID(i, j)] - r1 * F[ID(i + 1, j)])
                   + cc0[k] * (r1 * u[ID(i + 1, j)] - r0 * u[ID(i, j)]));
        double cl = CLAMP(cn, 0., 1.);
        ml += (cn - cl) * r;
        u_s[k] = cl;
      }
    mass_lost += ml * DX * DX;
  } else {
#pragma omp parallel for
    for (int i = 0; i < NR; i++)
      for (int j = 0; j <= NZ; j++) {
        double wf = w[ID(i, j)], s = wf * dt / DX;
        int donor = s > 0 ? ID(i, j - 1) : ID(i, j);
        F[ID(i, j)] = (wf == 0.) ? 0. : vof_flux_frac(donor, s, 1) * wf;
      }
    double ml = 0.0;
#pragma omp parallel for reduction(+ : ml)
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = ID(i, j);
        double cn = c[k] + dt / DX * ((F[ID(i, j)] - F[ID(i, j + 1)])
                   + cc0[k] * (w[ID(i, j + 1)] - w[ID(i, j)]));
        double cl = CLAMP(cn, 0., 1.);
        ml += (cn - cl) * rc(i);
        u_s[k] = cl;
      }
    mass_lost += ml * DX * DX;
  }
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double v = u_s[ID(i, j)];
      if (v < 1e-10) v = 0.; else if (v > 1. - 1e-10) v = 1.;
      c[ID(i, j)] = v;
    }
}

static void vof_advect(double dt, long step) {
  double *cc0 = w1; /* 作業配列を一時利用（移流前の c > 0.5 の指標）*/
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) cc0[ID(i, j)] = c[ID(i, j)] > 0.5 ? 1. : 0.;
  if (step % 2 == 0) { vof_sweep(0, dt, cc0); vof_sweep(1, dt, cc0); }
  else { vof_sweep(1, dt, cc0); vof_sweep(0, dt, cc0); }
  bc_c(c);
}

/* ------------------------------------------------------------------ 物性 */
static inline double rho_of(double cf) { return P.rho2 + (P.rho1 - P.rho2) * cf; }
static inline double mu_of(double cf) { return P.mu2 + (P.mu1 - P.mu2) * cf; }

static void update_props(void) {
  bc_c(c);
#pragma omp parallel for
  for (int i = -NG; i < NR + NG; i++)
    for (int j = -NG; j < NZ + NG; j++) {
      double cf = CLAMP(c[ID(i, j)], 0., 1.);
      rho_c[ID(i, j)] = rho_of(cf);
      mu_c[ID(i, j)] = mu_of(cf);
    }
}
static inline double rho_u(int i, int j) { return 0.5 * (rho_c[ID(i - 1, j)] + rho_c[ID(i, j)]); }
static inline double rho_w(int i, int j) { return 0.5 * (rho_c[ID(i, j - 1)] + rho_c[ID(i, j)]); }

/* ------------------------------------------------------------------ 曲率（高さ関数） */
static double hf_column_z(int i, int j, int below, int *ok) {
  /* 列 i で行 j-3..j+3 の和から界面の z 位置（中心セル中心からの Δ 単位）*/
  double H = 0.;
  for (int m = -3; m <= 3; m++) H += c[ID(i, j + m)];
  double cb = c[ID(i, j - 3)], ct = c[ID(i, j + 3)];
  if (below) { if (!(cb > 0.99 && ct < 0.01)) *ok = 0; return H - 3.5; }
  else { if (!(ct > 0.99 && cb < 0.01)) *ok = 0; return 3.5 - H; }
}
static double hf_column_r(int i, int j, int inside, int *ok) {
  double H = 0.;
  for (int m = -3; m <= 3; m++) H += c[ID(i + m, j)];
  double cl = c[ID(i - 3, j)], cr = c[ID(i + 3, j)];
  if (inside) { if (!(cl > 0.99 && cr < 0.01)) *ok = 0; return H - 3.5; }
  else { if (!(cr > 0.99 && cl < 0.01)) *ok = 0; return 3.5 - H; }
}

static void curvature(void) {
  bc_c(c);
  /* 平滑化 c（フォールバック CSF 用）: [1 2 1]^2/16 を 2 回 */
  for (int pass = 0; pass < 2; pass++) {
    const double *src = pass == 0 ? c : csm;
    double *dst = pass == 0 ? csm : csm2;
#pragma omp parallel for
    for (int i = -NG + 1; i < NR + NG - 1; i++)
      for (int j = -NG + 1; j < NZ + NG - 1; j++) {
        double s = 0.;
        for (int a = -1; a <= 1; a++)
          for (int b = -1; b <= 1; b++)
            s += (2 - abs(a)) * (2 - abs(b)) * src[ID(i + a, j + b)];
        dst[ID(i, j)] = s / 16.;
      }
  }
  double *cs = csm2;
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      int k = ID(i, j);
      kap[k] = 0.; hf_ok[k] = 0;
      double cc = c[k];
      if (cc <= 1e-6 || cc >= 1. - 1e-6) continue;
      double gr = (c[ID(i + 1, j + 1)] + 2 * c[ID(i + 1, j)] + c[ID(i + 1, j - 1)])
                - (c[ID(i - 1, j + 1)] + 2 * c[ID(i - 1, j)] + c[ID(i - 1, j - 1)]);
      double gz = (c[ID(i + 1, j + 1)] + 2 * c[ID(i, j + 1)] + c[ID(i - 1, j + 1)])
                - (c[ID(i + 1, j - 1)] + 2 * c[ID(i, j - 1)] + c[ID(i - 1, j - 1)]);
      int ok = 1;
      double kk = 0.;
      if (fabs(gz) >= fabs(gr)) {
        int below = gz < 0;
        double h0 = hf_column_z(i - 1, j, below, &ok);
        double h1 = hf_column_z(i, j, below, &ok);
        double h2 = hf_column_z(i + 1, j, below, &ok);
        if (ok && fabs(h1) <= 1.0) {
          double hp = 0.5 * (h2 - h0), hpp = h2 - 2 * h1 + h0;
          double q = 1. + hp * hp;
          double kp = hpp / (DX * pow(q, 1.5)) + hp / (rc(i) * sqrt(q));
          kk = below ? -kp : kp;
        } else ok = 0;
      } else {
        int inside = gr < 0;
        double h0 = hf_column_r(i, j - 1, inside, &ok);
        double h1 = hf_column_r(i, j, inside, &ok);
        double h2 = hf_column_r(i, j + 1, inside, &ok);
        double rint = rc(i) + h1 * DX;
        if (ok && fabs(h1) <= 1.0 && rint > 0.5 * DX) {
          double hp = 0.5 * (h2 - h0), hpp = h2 - 2 * h1 + h0;
          double q = 1. + hp * hp;
          double kp = 1. / (rint * sqrt(q)) - hpp / (DX * pow(q, 1.5));
          kk = inside ? kp : -kp;
        } else ok = 0;
      }
      if (ok) { kap[k] = CLAMP(kk, -2. / DX, 2. / DX); hf_ok[k] = 1; }
    }
  /* フォールバック: 近傍の有効な高さ関数曲率の平均、なければ平滑化 CSF */
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      int k = ID(i, j);
      double cc = c[k];
      if (hf_ok[k] || cc <= 1e-6 || cc >= 1. - 1e-6) continue;
      double s = 0.; int n = 0;
      for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++) {
          int ii = i + a, jj = j + b;
          if (ii < 0 || ii >= NR || jj < 0 || jj >= NZ) continue;
          if (hf_ok[ID(ii, jj)]) { s += kap[ID(ii, jj)]; n++; }
        }
      if (n > 0) { kap[k] = s / n; continue; }
      /* CSF: κ = -∇·(∇c̃/|∇c̃|)（軸対称）。面で法線を作って中心で発散 */
      double nrp, nrm, nzp, nzm;
      {
        double gr, gz, g;
        gr = (cs[ID(i + 1, j)] - cs[ID(i, j)]);
        gz = 0.25 * (cs[ID(i + 1, j + 1)] + cs[ID(i, j + 1)] - cs[ID(i + 1, j - 1)] - cs[ID(i, j - 1)]);
        g = sqrt(gr * gr + gz * gz) + 1e-12; nrp = gr / g;
        gr = (cs[ID(i, j)] - cs[ID(i - 1, j)]);
        gz = 0.25 * (cs[ID(i, j + 1)] + cs[ID(i - 1, j + 1)] - cs[ID(i, j - 1)] - cs[ID(i - 1, j - 1)]);
        g = sqrt(gr * gr + gz * gz) + 1e-12; nrm = gr / g;
        gz = (cs[ID(i, j + 1)] - cs[ID(i, j)]);
        gr = 0.25 * (cs[ID(i + 1, j + 1)] + cs[ID(i + 1, j)] - cs[ID(i - 1, j + 1)] - cs[ID(i - 1, j)]);
        g = sqrt(gr * gr + gz * gz) + 1e-12; nzp = gz / g;
        gz = (cs[ID(i, j)] - cs[ID(i, j - 1)]);
        gr = 0.25 * (cs[ID(i + 1, j)] + cs[ID(i + 1, j - 1)] - cs[ID(i - 1, j)] - cs[ID(i - 1, j - 1)]);
        g = sqrt(gr * gr + gz * gz) + 1e-12; nzm = gz / g;
      }
      double div = (rf(i + 1) * nrp - rf(i) * nrm) / (rc(i) * DX) + (nzp - nzm) / DX;
      kap[k] = CLAMP(-div, -2. / DX, 2. / DX);
    }
  bc_c(kap);
}

static inline double kappa_face(int k1, int k2) {
  double c1 = c[k1], c2 = c[k2];
  int i1 = c1 > 1e-6 && c1 < 1. - 1e-6, i2 = c2 > 1e-6 && c2 < 1. - 1e-6;
  if (i1 && i2) return 0.5 * (kap[k1] + kap[k2]);
  if (i1) return kap[k1];
  if (i2) return kap[k2];
  return 0.;
}

/* ------------------------------------------------------------------ 移流（運動量） */
static inline double mc_lim(double a, double b) {
  if (a * b <= 0.) return 0.;
  double s = a > 0 ? 1. : -1.;
  return s * MIN(MIN(2. * fabs(a), 2. * fabs(b)), 0.5 * fabs(a + b));
}
/* a ∂φ/∂x の風上 2 次（MC 制限）: f[-2..2] は φ の 5 点 */
static inline double upwind_deriv(double a, double fm2, double fm1, double f0, double fp1, double fp2) {
  double ph, mh;
  if (a > 0) {
    ph = f0 + 0.5 * mc_lim(f0 - fm1, fp1 - f0);
    mh = fm1 + 0.5 * mc_lim(fm1 - fm2, f0 - fm1);
  } else {
    ph = fp1 - 0.5 * mc_lim(fp1 - f0, fp2 - fp1);
    mh = f0 - 0.5 * mc_lim(f0 - fm1, fp1 - f0);
  }
  return a * (ph - mh) / DX;
}

static void advection_rhs(const double *U, const double *W, double *dU, double *dW) {
#pragma omp parallel for
  for (int i = 1; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double ua = U[ID(i, j)];
      double wa = 0.25 * (W[ID(i - 1, j)] + W[ID(i, j)] + W[ID(i - 1, j + 1)] + W[ID(i, j + 1)]);
      double ar = upwind_deriv(ua, U[ID(i - 2, j)], U[ID(i - 1, j)], U[ID(i, j)], U[ID(i + 1, j)], U[ID(i + 2, j)]);
      double az = upwind_deriv(wa, U[ID(i, j - 2)], U[ID(i, j - 1)], U[ID(i, j)], U[ID(i, j + 1)], U[ID(i, j + 2)]);
      dU[ID(i, j)] = -(ar + az);
    }
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 1; j <= NZ; j++) {
      double wa = W[ID(i, j)];
      double ua = 0.25 * (U[ID(i, j - 1)] + U[ID(i + 1, j - 1)] + U[ID(i, j)] + U[ID(i + 1, j)]);
      double ar = upwind_deriv(ua, W[ID(i - 2, j)], W[ID(i - 1, j)], W[ID(i, j)], W[ID(i + 1, j)], W[ID(i + 2, j)]);
      double az = upwind_deriv(wa, W[ID(i, j - 2)], W[ID(i, j - 1)], W[ID(i, j)], W[ID(i, j + 1)], W[ID(i, j + 2)]);
      dW[ID(i, j)] = -(ar + az);
    }
}

/* ------------------------------------------------------------------ 粘性 */
static inline double mu_corner(int i, int j) {
  /* r 面 i と z 面 j の交点（セル (i-1,j-1),(i,j-1),(i-1,j),(i,j) の角）*/
  return 0.25 * (mu_c[ID(i - 1, j - 1)] + mu_c[ID(i, j - 1)] + mu_c[ID(i - 1, j)] + mu_c[ID(i, j)]);
}
static inline double tau_rz(const double *U, const double *W, int i, int j) {
  return mu_corner(i, j) * ((U[ID(i, j)] - U[ID(i, j - 1)]) / DX + (W[ID(i, j)] - W[ID(i - 1, j)]) / DX);
}

static void viscous_add(const double *U, const double *W, double *Uo, double *Wo, double dt) {
#pragma omp parallel for
  for (int i = 1; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double r = rf(i);
      double srr_p = 2. * mu_c[ID(i, j)] * (U[ID(i + 1, j)] - U[ID(i, j)]) / DX * rc(i);
      double srr_m = 2. * mu_c[ID(i - 1, j)] * (U[ID(i, j)] - U[ID(i - 1, j)]) / DX * rc(i - 1);
      double muf = 0.5 * (mu_c[ID(i - 1, j)] + mu_c[ID(i, j)]);
      double f = (srr_p - srr_m) / (r * DX) - 2. * muf * U[ID(i, j)] / (r * r)
               + (tau_rz(U, W, i, j + 1) - tau_rz(U, W, i, j)) / DX;
      Uo[ID(i, j)] += dt * f / rho_u(i, j);
    }
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 1; j <= NZ; j++) {
      if (j == NZ && is_inflow_face(i)) continue;
      double r = rc(i);
      double szz_p = 2. * mu_c[ID(i, j)] * (W[ID(i, j + 1)] - W[ID(i, j)]) / DX;
      double szz_m = 2. * mu_c[ID(i, j - 1)] * (W[ID(i, j)] - W[ID(i, j - 1)]) / DX;
      double f = (rf(i + 1) * tau_rz(U, W, i + 1, j) - rf(i) * tau_rz(U, W, i, j)) / (r * DX)
               + (szz_p - szz_m) / DX;
      Wo[ID(i, j)] += dt * f / rho_w(i, j);
    }
}

/* ------------------------------------------------------------------ マルチグリッド */
typedef struct {
  int nr, nz, ni, nj;
  double *ax, *az, *dg, *di, *x, *b, *r;
} Level;
static Level *LV;
static int NLEV;
#define LID(L, i, j) (((i) + 1) * (L)->nj + (j) + 1)

static void mg_alloc(void) {
  NLEV = 1;
  int a = NR, b = NZ;
  while (a % 2 == 0 && b % 2 == 0 && a > 4 && b > 4) { a /= 2; b /= 2; NLEV++; }
  LV = calloc(NLEV, sizeof(Level));
  for (int l = 0; l < NLEV; l++) {
    Level *L = &LV[l];
    L->nr = NR >> l; L->nz = NZ >> l; L->ni = L->nr + 3; L->nj = L->nz + 3;
    size_t n = (size_t)L->ni * L->nj;
    L->ax = calloc(n, sizeof(double)); L->az = calloc(n, sizeof(double));
    L->dg = calloc(n, sizeof(double)); L->di = calloc(n, sizeof(double)); L->x = calloc(n, sizeof(double));
    L->b = calloc(n, sizeof(double)); L->r = calloc(n, sizeof(double));
  }
}

/* 細かいレベルの係数: ax(i,j) = r_f β_u, az(i,j) = r_c β_w, dg = 上端 Dirichlet の寄与 */
static void mg_setup(void) {
  Level *L = &LV[0];
#pragma omp parallel for
  for (int i = 0; i <= NR; i++)
    for (int j = 0; j <= NZ; j++) {
      int k = LID(L, i, j);
      L->ax[k] = (i > 0 && i < NR && j < NZ) ? rf(i) / rho_u(i, j) : 0.;
      L->az[k] = (j > 0 && j < NZ && i < NR) ? rc(i) / rho_w(i, j) : 0.;
      L->dg[k] = 0.;
      if (j == NZ - 1 && i < NR && !is_inflow_face(i)) L->dg[k] = 2. * rc(i) / rho_c[ID(i, NZ - 1)];
    }
  for (int l = 1; l < NLEV; l++) {
    Level *F = &LV[l - 1], *C = &LV[l];
#pragma omp parallel for
    for (int I = 0; I <= C->nr; I++)
      for (int J = 0; J <= C->nz; J++) {
        int k = LID(C, I, J);
        C->ax[k] = (I <= C->nr && J < C->nz) ? F->ax[LID(F, 2 * I, 2 * J)] + F->ax[LID(F, 2 * I, 2 * J + 1)] : 0.;
        C->az[k] = (J <= C->nz && I < C->nr) ? F->az[LID(F, 2 * I, 2 * J)] + F->az[LID(F, 2 * I + 1, 2 * J)] : 0.;
        C->dg[k] = (I < C->nr && J < C->nz)
                     ? F->dg[LID(F, 2 * I, 2 * J)] + F->dg[LID(F, 2 * I + 1, 2 * J)]
                       + F->dg[LID(F, 2 * I, 2 * J + 1)] + F->dg[LID(F, 2 * I + 1, 2 * J + 1)]
                     : 0.;
      }
  }
}

static void mg_diag(void) {
  for (int l = 0; l < NLEV; l++) {
    Level *L = &LV[l];
#pragma omp parallel for if (L->nr * L->nz > 8192)
    for (int i = 0; i < L->nr; i++)
      for (int j = 0; j < L->nz; j++) {
        int k = LID(L, i, j);
        double d = L->ax[k] + L->ax[LID(L, i + 1, j)] + L->az[k] + L->az[LID(L, i, j + 1)] + L->dg[k];
        L->di[k] = d > 0 ? 1. / d : 0.;
      }
  }
}

static inline double apply_A(Level *L, const double *x, int i, int j) {
  int k = LID(L, i, j);
  double xc = x[k];
  double s = L->ax[k] * (xc - x[LID(L, i - 1, j)]) + L->ax[LID(L, i + 1, j)] * (xc - x[LID(L, i + 1, j)])
           + L->az[k] * (xc - x[LID(L, i, j - 1)]) + L->az[LID(L, i, j + 1)] * (xc - x[LID(L, i, j + 1)])
           + L->dg[k] * xc;
  return s;
}

static void smooth(Level *L, int color) {
#pragma omp parallel for if (L->nr * L->nz > 8192)
  for (int i = 0; i < L->nr; i++)
    for (int j = (i + color) & 1; j < L->nz; j += 2) {
      int k = LID(L, i, j), nj = L->nj;
      L->x[k] = (L->b[k] + L->ax[k] * L->x[k - nj] + L->ax[k + nj] * L->x[k + nj]
                 + L->az[k] * L->x[k - 1] + L->az[k + 1] * L->x[k + 1]) * L->di[k];
    }
}

static void vcycle(int l) {
  Level *L = &LV[l];
  if (l == NLEV - 1) {
    for (int it = 0; it < 60; it++) { smooth(L, 0); smooth(L, 1); smooth(L, 1); smooth(L, 0); }
    return;
  }
  for (int s = 0; s < P.mg_pre; s++) { smooth(L, 0); smooth(L, 1); }
#pragma omp parallel for if (L->nr * L->nz > 8192)
  for (int i = 0; i < L->nr; i++)
    for (int j = 0; j < L->nz; j++) L->r[LID(L, i, j)] = L->b[LID(L, i, j)] - apply_A(L, L->x, i, j);
  Level *C = &LV[l + 1];
#pragma omp parallel for if (C->nr * C->nz > 8192)
  for (int I = 0; I < C->nr; I++)
    for (int J = 0; J < C->nz; J++) {
      int k = LID(C, I, J);
      C->b[k] = L->r[LID(L, 2 * I, 2 * J)] + L->r[LID(L, 2 * I + 1, 2 * J)]
              + L->r[LID(L, 2 * I, 2 * J + 1)] + L->r[LID(L, 2 * I + 1, 2 * J + 1)];
      C->x[k] = 0.;
    }
  vcycle(l + 1);
#pragma omp parallel for if (L->nr * L->nz > 8192)
  for (int i = 0; i < L->nr; i++)
    for (int j = 0; j < L->nz; j++) L->x[LID(L, i, j)] += C->x[LID(C, i / 2, j / 2)];
  for (int s = 0; s < P.mg_post; s++) { smooth(L, 1); smooth(L, 0); }
}

/* Flexible PCG: A x = b（細かいレベル）。返り値は反復回数 */
static double *cg_r, *cg_z, *cg_zo, *cg_d, *cg_q;
static int pcg_solve(double *x, const double *b, double tol_abs, double *res_out) {
  Level *L = &LV[0];
  size_t n = (size_t)L->ni * L->nj;
  if (!cg_r) {
    cg_r = calloc(n, sizeof(double)); cg_z = calloc(n, sizeof(double)); cg_zo = calloc(n, sizeof(double));
    cg_d = calloc(n, sizeof(double)); cg_q = calloc(n, sizeof(double));
  }
  double rmax = 0.;
#pragma omp parallel for reduction(max : rmax)
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      int k = LID(L, i, j);
      cg_r[k] = b[k] - apply_A(L, x, i, j);
      double v = fabs(cg_r[k]) / rc(i);
      if (v > rmax) rmax = v;
    }
  if (rmax < tol_abs) { *res_out = rmax; return 0; }
  double rz_old = 0.;
  int it;
  for (it = 1; it <= 400; it++) {
    /* z = M^-1 r */
    memcpy(L->b, cg_r, n * sizeof(double));
    memset(L->x, 0, n * sizeof(double));
    vcycle(0);
    double rz = 0., rzo = 0.;
#pragma omp parallel for reduction(+ : rz, rzo)
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = LID(L, i, j);
        rz += cg_r[k] * L->x[k];
        rzo += cg_r[k] * cg_zo[k];
      }
    double beta = (it == 1) ? 0. : (rz - rzo) / rz_old;  /* Polak–Ribière（flexible）*/
    if (beta < 0) beta = 0.;
#pragma omp parallel for
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = LID(L, i, j);
        cg_d[k] = L->x[k] + beta * cg_d[k];
        cg_zo[k] = L->x[k];
      }
    rz_old = rz;
    double dq = 0.;
#pragma omp parallel for reduction(+ : dq)
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = LID(L, i, j);
        cg_q[k] = apply_A(L, cg_d, i, j);
        dq += cg_d[k] * cg_q[k];
      }
    double alpha = rz / dq;
    rmax = 0.;
#pragma omp parallel for reduction(max : rmax)
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k = LID(L, i, j);
        x[k] += alpha * cg_d[k];
        cg_r[k] -= alpha * cg_q[k];
        double v = fabs(cg_r[k]) / rc(i);
        if (v > rmax) rmax = v;
      }
    if (rmax < tol_abs) break;
  }
  *res_out = rmax;
  return it;
}

/* ------------------------------------------------------------------ 射影 */
static double *pres_mg; /* MG レイアウトの圧力（前ステップを初期値に使う）*/
static int last_iters = 0;
static double last_res = 0.;

static void project(double dt) {
  bc_vel(u, w);
  mg_setup();
  mg_diag();
  Level *L = &LV[0];
  double *b = calloc((size_t)L->ni * L->nj, sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double dv = rf(i + 1) * u[ID(i + 1, j)] - rf(i) * u[ID(i, j)] + rc(i) * (w[ID(i, j + 1)] - w[ID(i, j)]);
      b[LID(L, i, j)] = -DX * dv / dt;
    }
  /* 許容残差: 発散 × dt が tol 以下（|r|/r_c = Δ^2 |div|/dt）*/
  double tol_abs = P.tol * DX * DX / (dt * dt);
  last_iters = pcg_solve(pres_mg, b, tol_abs, &last_res);
  last_res = last_res * dt * dt / (DX * DX);
  free(b);
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) p[ID(i, j)] = pres_mg[LID(L, i, j)];
  /* 速度修正 */
#pragma omp parallel for
  for (int i = 1; i < NR; i++)
    for (int j = 0; j < NZ; j++)
      u[ID(i, j)] -= dt / rho_u(i, j) * (p[ID(i, j)] - p[ID(i - 1, j)]) / DX;
#pragma omp parallel for
  for (int i = 0; i < NR; i++) {
    for (int j = 1; j < NZ; j++)
      w[ID(i, j)] -= dt / rho_w(i, j) * (p[ID(i, j)] - p[ID(i, j - 1)]) / DX;
    if (!is_inflow_face(i))
      w[ID(i, NZ)] -= dt / rho_c[ID(i, NZ - 1)] * (0. - p[ID(i, NZ - 1)]) / (0.5 * DX);
  }
  bc_vel(u, w);
}

/* ------------------------------------------------------------------ 初期条件 */
/* 高さ z での初期水柱の速度と半径（上端で jet_velocity, jet_radius）*/
static double jet_v_at(double z) {
  double Lz = NZ * DX;
  if (!P.jet_profile) return P.jet_velocity;
  return sqrt(SQ(P.jet_velocity) + 2. * P.grav * MAX(Lz - z, 0.));
}
static double jet_a_at(double z) {
  return P.jet_radius * sqrt(P.jet_velocity / jet_v_at(z));
}
/* 点 (r,z) が液体か */
static int liquid_at(double r, double z) {
  int liq = 0;
  if (P.pool_depth > 0 && z < P.pool_depth) liq = 1;
  if (P.jet_radius > 0 && P.jet_init) {
    double a = jet_a_at(z), R = P.jet_head_radius > 0 ? P.jet_head_radius : jet_a_at(P.jet_tip_z);
    double zc = P.jet_tip_z + R; /* 先頭の球の中心 */
    if (z >= zc && r < a) liq = 1;
    if (SQ(r) + SQ(z - zc) < R * R) liq = 1;
  }
  if (P.drop_radius > 0) {
    double dz = z - P.drop_z, rr = sqrt(r * r + dz * dz), ct = rr > 0 ? dz / rr : 1.;
    double Rd = P.drop_radius * (1. + P.drop_p2 * 0.5 * (3. * ct * ct - 1.));
    if (rr < Rd) liq = 1;
  }
  if (P.bubble_radius > 0 && SQ(r) + SQ(z - P.bubble_z) < SQ(P.bubble_radius)) liq = 0;
  if (P.cavity_radius > 0) {
    /* 水面直下に半径 cavity_radius の球状のくぼみ（中心は pool_depth - cavity_depth）*/
    double zc = P.pool_depth - P.cavity_depth;
    if (SQ(r) + SQ(z - zc) < SQ(P.cavity_radius) && z < P.pool_depth + 1e-12) liq = 0;
  }
  return liq;
}
static double initial_w(double r, double z) {
  if (P.drop_radius > 0 && SQ(r) + SQ(z - P.drop_z) < SQ(P.drop_radius * 1.02)) return -P.drop_velocity;
  if (P.jet_radius > 0 && P.jet_init) {
    double a = jet_a_at(z), R = P.jet_head_radius > 0 ? P.jet_head_radius : jet_a_at(P.jet_tip_z);
    double zc = P.jet_tip_z + R;
    if (z >= zc && r < a * 1.05) return -jet_v_at(z);
    if (SQ(r) + SQ(z - zc) < SQ(R * 1.05)) return -jet_v_at(zc);
  }
  return 0.;
}

static void init_fields(void) {
  const int S = 12;
#pragma omp parallel for
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double sum = 0., wsum = 0.;
      for (int a = 0; a < S; a++)
        for (int b = 0; b < S; b++) {
          double r = (i + (a + 0.5) / S) * DX, z = (j + (b + 0.5) / S) * DX;
          sum += r * liquid_at(r, z);
          wsum += r;
        }
      c[ID(i, j)] = sum / wsum;
    }
#pragma omp parallel for
  for (int i = 0; i <= NR; i++)
    for (int j = 0; j <= NZ; j++) {
      w[ID(i, j)] = (i < NR) ? initial_w(rc(i), rf(j)) : 0.;
      u[ID(i, j)] = 0.;
    }
  bc_c(c);
  bc_vel(u, w);
}

/* ------------------------------------------------------------------ 出力 */
static void write_snapshot(int idx) {
  char fn[700];
  snprintf(fn, sizeof fn, "%s/snap_%05d.bin", P.out_dir, idx);
  FILE *f = fopen(fn, "wb");
  if (!f) { perror(fn); return; }
  int hdr[2] = {NR, NZ};
  double hd[2] = {DX, t_now};
  fwrite(hdr, sizeof(int), 2, f);
  fwrite(hd, sizeof(double), 2, f);
  float *buf = malloc(sizeof(float) * NR * NZ);
  /* c, セル中心 u, w, p の順 */
  for (int q = 0; q < 4; q++) {
    for (int i = 0; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        double v;
        if (q == 0) v = c[ID(i, j)];
        else if (q == 1) v = 0.5 * (u[ID(i, j)] + u[ID(i + 1, j)]);
        else if (q == 2) v = 0.5 * (w[ID(i, j)] + w[ID(i, j + 1)]);
        else v = p[ID(i, j)];
        buf[i * NZ + j] = (float)v;
      }
    fwrite(buf, sizeof(float), NR * NZ, f);
  }
  free(buf);
  fclose(f);
}

static double total_volume(void) {
  double s = 0.;
#pragma omp parallel for reduction(+ : s)
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) s += c[ID(i, j)] * rc(i);
  return 2. * M_PI * s * DX * DX;
}

static void stats(double *umax, double *ke) {
  double m = 0., e = 0.;
#pragma omp parallel for reduction(max : m) reduction(+ : e)
  for (int i = 0; i < NR; i++)
    for (int j = 0; j < NZ; j++) {
      double uu = 0.5 * (u[ID(i, j)] + u[ID(i + 1, j)]), ww = 0.5 * (w[ID(i, j)] + w[ID(i, j + 1)]);
      double v = sqrt(uu * uu + ww * ww);
      if (v > m) m = v;
      e += 0.5 * rho_c[ID(i, j)] * (uu * uu + ww * ww) * rc(i);
    }
  *umax = m; *ke = 2. * M_PI * e * DX * DX;
}

/* ------------------------------------------------------------------ メイン */
int main(int argc, char **argv) {
  params_default(&P);
  if (argc > 1) params_read(&P, argv[1]);
  NR = P.nr; NZ = P.nz; DX = P.Lr / NR;
  NI = NR + 2 * NG + 1; NJ = NZ + 2 * NG + 1;
  char cmd[700];
  snprintf(cmd, sizeof cmd, "mkdir -p %s", P.out_dir);
  if (system(cmd)) {}
  c = alloc_field(); u = alloc_field(); w = alloc_field(); p = alloc_field(); kap = alloc_field();
  mx_ = alloc_field(); my_ = alloc_field(); alph = alloc_field();
  u_s = alloc_field(); w_s = alloc_field(); u1 = alloc_field(); w1 = alloc_field();
  du = alloc_field(); dw = alloc_field(); rho_c = alloc_field(); mu_c = alloc_field();
  csm = alloc_field(); csm2 = alloc_field();
  hf_ok = calloc((size_t)NI * NJ, sizeof(int));
  mg_alloc();
  pres_mg = calloc((size_t)LV[0].ni * LV[0].nj, sizeof(double));

  init_fields();
  update_props();
  curvature();
  double V0 = total_volume();
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  fprintf(stderr, "grid %d x %d, dx = %.3g mm, Lz = %.3g mm, MG levels %d, threads %d\n",
          NR, NZ, DX * 1e3, NZ * DX * 1e3, NLEV, nthreads);

  char fn[700];
  snprintf(fn, sizeof fn, "%s/diag.csv", P.out_dir);
  FILE *fd = fopen(fn, "w");
  fprintf(fd, "step,t,dt,iters,res,volume,mass_err,umax,ke,wall\n");

  int snap = 0;
  double next_out = 0.;
  write_snapshot(snap++);
  next_out += P.out_dt;
  double dt_cap = sqrt((P.rho1 + P.rho2) * DX * DX * DX / (4. * M_PI * P.sigma));
  double nu_max = MAX(P.mu1 / P.rho1, P.mu2 / P.rho2);
  double dt_visc = 0.2 * DX * DX / nu_max;
  double wt0 = 0;
#ifdef _OPENMP
  wt0 = omp_get_wtime();
#endif
  long step;
  double dt = 0.;
  for (step = 0; step < P.max_steps && t_now < P.t_end - 1e-12; step++) {
    /* --- 時間刻み --- */
    double umax = 1e-6;
#pragma omp parallel for reduction(max : umax)
    for (int i = 0; i <= NR; i++)
      for (int j = 0; j <= NZ; j++) {
        double a = fabs(u[ID(i, j)]), b = fabs(w[ID(i, j)]);
        if (a > umax) umax = a;
        if (b > umax) umax = b;
      }
    double dtn = MIN(P.cfl * DX / umax, 0.8 * dt_cap);
    dtn = MIN(dtn, dt_visc);
    if (dt > 0) dtn = MIN(dtn, 1.2 * dt);
    if (t_now + dtn > next_out) dtn = MAX(next_out - t_now, 1e-12);
    if (t_now + dtn > P.t_end) dtn = P.t_end - t_now;
    dt = dtn;

    /* --- 1. 界面移流 --- */
    double tA = wclock();
    vof_advect(dt, step);
    double tB = wclock();
    update_props();
    curvature();
    double tC = wclock();

    /* --- 2. 運動量: 移流（SSP-RK2）--- */
    bc_vel(u, w);
    memcpy(u_s, u, sizeof(double) * NI * NJ);
    memcpy(w_s, w, sizeof(double) * NI * NJ);
    memset(du, 0, sizeof(double) * NI * NJ);
    memset(dw, 0, sizeof(double) * NI * NJ);
    advection_rhs(u, w, du, dw);
#pragma omp parallel for
    for (int k = 0; k < NI * NJ; k++) { u1[k] = u[k] + dt * du[k]; w1[k] = w[k] + dt * dw[k]; }
    bc_vel(u1, w1);
    advection_rhs(u1, w1, du, dw);
#pragma omp parallel for
    for (int k = 0; k < NI * NJ; k++) {
      u_s[k] = 0.5 * (u[k] + u1[k] + dt * du[k]);
      w_s[k] = 0.5 * (w[k] + w1[k] + dt * dw[k]);
    }
    /* 粘性（u^n で評価）・重力・表面張力 */
    viscous_add(u, w, u_s, w_s, dt);
#pragma omp parallel for
    for (int i = 1; i < NR; i++)
      for (int j = 0; j < NZ; j++) {
        int k1 = ID(i - 1, j), k2 = ID(i, j);
        double fst = P.sigma * kappa_face(k1, k2) * (c[k2] - c[k1]) / DX;
        u_s[ID(i, j)] += dt * fst / rho_u(i, j);
      }
#pragma omp parallel for
    for (int i = 0; i < NR; i++)
      for (int j = 1; j <= NZ; j++) {
        if (j == NZ && is_inflow_face(i)) continue;
        double fst = 0.;
        if (j < NZ) {
          int k1 = ID(i, j - 1), k2 = ID(i, j);
          fst = P.sigma * kappa_face(k1, k2) * (c[k2] - c[k1]) / DX;
        }
        double rho = (j < NZ) ? rho_w(i, j) : rho_c[ID(i, NZ - 1)];
        w_s[ID(i, j)] += dt * (fst / rho - P.grav);
      }
    memcpy(u, u_s, sizeof(double) * NI * NJ);
    memcpy(w, w_s, sizeof(double) * NI * NJ);

    /* --- 3. 圧力射影 --- */
    double tD = wclock();
    project(dt);
    double tE = wclock();
    T_vof += tB - tA; T_curv += tC - tB; T_mom += tD - tC; T_proj += tE - tD;
    t_now += dt;

    /* --- 出力 --- */
    if (t_now >= next_out - 1e-12 || t_now >= P.t_end - 1e-12) {
      double um, ke;
      stats(&um, &ke);
      double V = total_volume();
      double wt = 0;
#ifdef _OPENMP
      wt = omp_get_wtime() - wt0;
#else
      wt = 0.;
#endif
      fprintf(fd, "%ld,%.8e,%.4e,%d,%.3e,%.10e,%.3e,%.4e,%.6e,%.2f\n", step, t_now, dt, last_iters,
              last_res, V, mass_lost * 2 * M_PI / V0, um, ke, wt);
      fflush(fd);
      if (P.verbose)
        fprintf(stderr, "step %6ld t=%.5f dt=%.2e it=%3d res=%.1e V/V0-1=%+.2e umax=%.3f wall=%.0fs\n",
                step, t_now, dt, last_iters, last_res, V / V0 - 1., um, wt);
      write_snapshot(snap++);
      next_out += P.out_dt;
    }
  }
  fclose(fd);
  fprintf(stderr, "done: %ld steps  [time: vof %.2fs, props+curv %.2fs, momentum %.2fs, projection %.2fs]\n",
          step, T_vof, T_curv, T_mom, T_proj);
  return 0;
}
