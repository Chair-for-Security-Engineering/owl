#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef double   f64;
#include "params.c"

#define LEN(a)     (sizeof(a) / sizeof(*a))
#define LEN_LOGN   LEN(params_logsN)
#define LEN_SINGLE LEN(params_single)
#define LEN_FOLDED LEN(params_folded)
#define LEN_DOUBLE LEN(params_double)
#define LEN_SWITCH LEN(params_switch)
#define LEN_PARAMS (LEN_SINGLE + LEN_FOLDED + LEN_DOUBLE + LEN_SWITCH)

#define TABLE(lbl, cfg)                           \
"\\begin{table}\n"                                \
"    \\centering\\setlength{\\tabcolsep}{.5em}\n" \
"    \\caption{\\the\\owlevalcap}\n"              \
"    \\label{owl:tbl:" #lbl "}\n"                 \
"    \\begin{tabular}{" #cfg "}\n"                \
"    \\toprule\n"
#define TLONG(lbl, cfg)                           \
"\\begin{longtable}{" #cfg "}\n"                  \
"\\caption{\\the\\owlevalcap}\n"                  \
"\\label{owl:tbl:" #lbl "}\n"                     \
"\\\\\\toprule\n"
#define TABLE_MID                                 \
"    \\midrule\n"
#define TABLE_END                                 \
"    \\bottomrule\n"                              \
"    \\end{tabular}\n"                            \
"\\end{table}\n\n"
#define TLONG_END                                 \
"\\bottomrule\n"                                  \
"\\end{longtable}\n\n"

static int
cmp(void const *op1, void const *op2) {
	f64 t1 = *(f64 *)op1;
	f64 t2 = *(f64 *)op2;
	return (t1 > t2) - (t1 < t2);
}

void
printq1(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = *idx0;
	fputs(TABLE(omega, crrr)
	      "    \\multicolumn{2}{c}{Parameters} & \\multicolumn{2}{c}{Time (ms)} \\\\\n"
	      "    \\cmidrule(lr{1pt}){1-2}\\cmidrule(lr{1pt}){3-4} $\\log_2 N$ & \\multicolumn{1}{c}{$\\log_2 q$} & $\\omega = 1$ & $\\omega = 2$ \\\\\n"
	      TABLE_MID, f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		f64 maxQP = (f64)params_logsQP[idx];
		Params *p = params_folded + pidx;
		f64 *t = times + LEN_SINGLE + pidx;

		u32 L = (u32)p[0].L;
		f64 logQ = 0.0;
		for (u32 i = 0; i < L; ++i)
			logQ += log2((f64)p[0].mods[i]);
		f64 t0 = t[0] * 1000.0, t1 = t[1] * 1000.0;
		#define Q1_FMT "    %2d & %4ld (%5.2f\\,\\%%) & %6.1f & %6.1f \\\\\n"
		fprintf(f, Q1_FMT, logN, lround(logQ), logQ / maxQP * 100.0, t0, t1);

		pidx += 2;
	}
	fputs(TABLE_END, f);
	*idx0 = pidx;
}

void
printq2(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = *idx0;
	fputs(TABLE(degree, crrcrrc)
	      "    \\multicolumn{4}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "    \\cmidrule(lr{1pt}){1-4}\\cmidrule(lr{1pt}){5-7} $\\log_2 N$ & \\multicolumn{1}{c}{$\\log_2 q$} & $\\omega$ & $\\omega'$ & $N$ & $2 N$ & Speedup \\\\\n"
	      TABLE_MID, f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		f64 maxQP = (f64)params_logsQP[idx];
		Params *p = params_folded + pidx;
		f64 *t = times + LEN_SINGLE + pidx;

		u32 L0 = (u32)p[0].L;
		int W0 = (int)p[0].W;
		f64 logQ0 = 0.0;
		for (u32 i = 0; i < L0; ++i)
			logQ0 += log2((f64)p[0].mods[i]);
		u32 L2 = (u32)p[2].L;
		int W2 = (int)p[2].W;
		f64 logQ2 = 0.0;
		for (u32 i = 0; i < L2; ++i)
			logQ2 += log2((f64)p[0].mods[i]);
		f64 t0 = t[0] * 1000.0, t1 = t[1] * 1000.0;
		f64 t2 = t[2] * 1000.0, t3 = t[3] * 1000.0;

		#define Q2_FMT "    %2d & %4ld (%5.2f\\,\\%%) & %2d & %1d & %6.1f & %6.1f & %4.2f \\\\\n"
		fprintf(f, Q2_FMT, logN, lround(logQ0), logQ0 / maxQP * 100.0, W0, 1, t0, t1, t0 / t1);
		fprintf(f, Q2_FMT, logN, lround(logQ2), logQ2 / maxQP * 100.0, W2, 1, t2, t3, t2 / t3);
		if (idx != LEN_LOGN - 1) fputs(TABLE_MID, f);

		pidx += 4;
	}
	fputs(TABLE_END, f);
	*idx0 = pidx;
}

void
printq3(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = *idx0;
	f64 xshift = -3.2;
	fputs("\\begin{figure}\n" \
	      "    \\centering\n" \
	      "    \\begin{tikzpicture}[>=latex]\n", f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		fprintf(f, "    \\node at (%.1f, %.1f) {$N = 2^{%d}$};\n", xshift, 2.5, (int)params_logsN[idx]);
		xshift += 3.0;
	}
	xshift = -4.5;
	u32 tidx = LEN_SINGLE + LEN_FOLDED;
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		Params *p = params_folded + pidx;
		f64 *tfld = times + LEN_SINGLE + pidx;
		f64 *tdbl = times + tidx;
		u32 L = (u32)p[0].L;

		#define Q3_SCOPE                                                   \
		"    \\begin{scope}[shift={(%.2f, %.2f)}, xscale=.02, yscale=2]\n" \
		"        \\draw[->] (0, 0) -- ++(100, 0) node[right] {$\\ell$};\n" \
		"        \\draw[->] (0, 0) -- ++(0, 1) node[left] {t};\n"          \
		"        \\draw[green]\n"
		#define Q3_SCOPE_MID                                               \
		"        \\draw[blue]\n"
		#define Q3_SCOPE_END                                               \
		"    \\end{scope}\n"
		#define Q3_SCOPE_ROW                                               \
		"            (%5.1f, %7.5f) node[dot=.1] {} %s\n"

		f64 tmax = 0.0;
		for (u32 i = 0; i < L; ++i) {
			tmax = tmax > tfld[i] ? tmax : tfld[i];
			tmax = tmax > tdbl[i] ? tmax : tdbl[i];
			tmax = tmax > tfld[L + i] ? tmax : tfld[L + i];
			tmax = tmax > tdbl[L + i] ? tmax : tdbl[L + i];
		}
		fprintf(f, Q3_SCOPE, xshift, 0.0);
		for (u32 i = 0; i < L; ++i) {
			f64 t = tdbl[i] / tmax;
			char const *s = i == L - 1 ? ";" : "--";
			fprintf(f, Q3_SCOPE_ROW, (f64)(L - i) / (f64)L * 100.0, t, s);
		}
		fputs(Q3_SCOPE_MID, f);
		for (u32 i = 0; i < L; ++i) {
			f64 t = tfld[i] / tmax;
			char const *s = i == L - 1 ? ";" : "--";
			fprintf(f, Q3_SCOPE_ROW, (f64)(L - i) / (f64)L * 100.0, t, s);
		}
		fprintf(f, Q3_SCOPE_END Q3_SCOPE, xshift, -3.0);
		for (u32 i = 0; i < L; ++i) {
			f64 t = tdbl[L + i] / tmax;
			char const *s = i == L - 1 ? ";" : "--";
			fprintf(f, Q3_SCOPE_ROW, (f64)(L - i) / (f64)L * 100.0, t, s);
		}
		fputs(Q3_SCOPE_MID, f);
		for (u32 i = 0; i < L; ++i) {
			f64 t = tfld[L + i] / tmax;
			char const *s = i == L - 1 ? ";" : "--";
			fprintf(f, Q3_SCOPE_ROW, (f64)(L - i) / (f64)L * 100.0, t, s);
		}
		fputs(Q3_SCOPE_END, f);

		pidx += 2 * L, tidx += 2 * L;
		xshift += 3.0;
	}
	fputs("    \\end{tikzpicture}\n" \
	      "    \\caption{\\the\\owlevalcap}\n" \
	      "    \\label{owl:fig:decomposition}\n" \
	      "\\end{figure}\n\n", f);

	pidx = 0;
	fputs(TLONG(decomposition, lcrrrrr)
	      "Type & $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & \\multicolumn{1}{c}{$k$} & \\multicolumn{1}{c}{$\\omega$} & \\multicolumn{1}{c}{$\\tilde\\omega$} & \\multicolumn{1}{c}{Time (ms)} \\\\\n"
	      "\\midrule\\endfirsthead\n"
	      "Type & $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & \\multicolumn{1}{c}{$k$} & \\multicolumn{1}{c}{$\\omega$} & \\multicolumn{1}{c}{$\\tilde\\omega$} & \\multicolumn{1}{c}{Time (ms)} \\\\\n"
	      "\\midrule\\endhead\n", f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		Params *pfld = params_folded + *idx0 + pidx;
		Params *pdbl = params_double + pidx;
		f64 *tfld = times + LEN_SINGLE + *idx0 + pidx;
		f64 *tdbl = times + LEN_SINGLE + LEN_FOLDED + pidx;

		u32 L = (u32)pfld[0].L;
		for (u32 i = 0; i < L; ++i) {
			#define Q3_FMT "%s & %2d & %2d & %2d & %2d & %2d & %6.1f \\\\\n"
			fprintf(f, Q3_FMT, "single", logN, (int)(L - i), 1, (int)(L - i), 0, tfld[i] * 1000.0);
			fprintf(f, Q3_FMT, "single", logN, (int)(L - i), (int)pfld[L + i].K, (int)pfld[L + i].W, 0, tfld[L + i] * 1000.0);
			fprintf(f, Q3_FMT, "double", logN, (int)(L - i), 1, (int)(L - i), (int)pdbl[i].W2, tdbl[i] * 1000.0);
			fprintf(f, Q3_FMT, "double", logN, (int)(L - i), (int)pdbl[L + i].K, (int)pdbl[L + i].W, (int)pdbl[L + i].W2, tdbl[L + i] * 1000.0);
		}

		pidx += 2 * L, tidx += 2 * L;
	}
	fputs(TLONG_END, f);
	*idx0 += pidx;
}

void
printq4(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = *idx0, cnt = 0;
	fputs(TLONG(primesize, crrrrrr)
	      "\\multicolumn{4}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "\\cmidrule(lr{1pt}){1-4}\\cmidrule(lr{1pt}){5-7} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell_{36}$} & \\multicolumn{1}{c}{$\\ell_{54}$} & $\\omega$ & \\texttt{small} & \\texttt{large} & Speedup \\\\\n"
	      "\\midrule\\endfirsthead\n"
	      "\\cmidrule(lr{1pt}){1-4}\\cmidrule(lr{1pt}){5-7} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell_{36}$} & \\multicolumn{1}{c}{$\\ell_{54}$} & $\\omega$ & \\texttt{small} & \\texttt{large} & Speedup \\\\\n"
	      "\\midrule\\endhead\n", f);
	f64 avg = 0.0, absL[LEN_LOGN];
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		Params *p = params_folded + pidx;
		f64 *t = times + LEN_SINGLE + pidx;

		u32 L = (u32)p[1].L;
		for (u32 i = 0; i < L / 2; ++i) {
			int L36 = (int)p[2 * i].L;
			int L54 = (int)p[2 * i + 1].L;
			int W = (int)p[2 * i + 1].W;
			f64 t0 = t[2 * i] * 1000.0, t1 = t[2 * i + 1] * 1000.0;
			#define Q4_FMT "%2d & %2d & %2d & %2d & %6.1f & %6.1f & %4.2f \\\\\n"
			fprintf(f, Q4_FMT, logN, L36, L54, W, t0, t1, t0 / t1);
			avg += t0 / t1, ++cnt;
			if (i == 0)
				absL[idx] = t1 - t0;
		}
		pidx += L;
	}
	fputs(TLONG_END, f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx)
		fprintf(f, "[+] absolute: %lf\n", absL[idx]);
	fprintf(f, "[+] average:  %lf\n\n", avg / (f64)cnt);
}

void
printq5(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = *idx0;
	u32 tidx = LEN_SINGLE + LEN_FOLDED + LEN_DOUBLE;
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		Params *p = params_folded + pidx;
		f64 *tfld = times + LEN_SINGLE + pidx;
		f64 *tsw = times + tidx;

		u32 L = (u32)p[1].L;
		f64 diff[L / 2];
		for (u32 i = 0; i < L / 2; ++i)
			diff[i] = (tsw[i] - tfld[2 * i + 1]) * 1000.0;
		qsort(diff, L / 2, 8, cmp);
		fprintf(f, "[+] overhead 2^%2d: %lf\n", logN, diff[L / 4]);
		pidx += L, tidx += L / 2;
	}

	pidx = *idx0;
	fputs("\n" TLONG(replacing, crrrrr)
	      "\\multicolumn{3}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "\\cmidrule(lr{1pt}){1-3}\\cmidrule(lr{1pt}){4-6} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & $\\omega$ & \\texttt{folded} & \\texttt{switch} & Overhead \\\\\n"
	      "\\midrule\\endfirsthead\n"
	      "\\multicolumn{3}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "\\cmidrule(lr{1pt}){1-3}\\cmidrule(lr{1pt}){4-6} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & $\\omega$ & \\texttt{folded} & \\texttt{switch} & Overhead \\\\\n"
	      "\\midrule\\endhead\n", f);
	tidx = LEN_SINGLE + LEN_FOLDED + LEN_DOUBLE;
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		Params *p = params_folded + pidx;
		f64 *tfld = times + LEN_SINGLE + pidx;
		f64 *tsw = times + tidx;

		u32 L = (u32)p[1].L;
		for (u32 i = 0; i < L / 2; ++i) {
			int W = (int)p[2 * i + 1].W;
			f64 t0 = tfld[2 * i + 1] * 1000.0, t1 = tsw[i] * 1000.0;
			#define Q5_FMT "%2d & %2d & %2d & %6.1f & %6.1f & %5.1f \\\\\n"
			fprintf(f, Q5_FMT, logN, (int)(L - 2 * i), W, t0, t1, t1 - t0);
		}
		pidx += L, tidx += L / 2;
	}
	fputs(TLONG_END, f);
	*idx0 = pidx;
}

void
printq6(FILE *f, f64 *times, u32 *idx0)
{
	u32 pidx = 0;
	fputs(TABLE(folding, cccccc)
	      "    $\\log_2 N$ & \\multicolumn{5}{c}{$\\omega$} \\\\\n" \
	      "    \\cmidrule(lr{1pt}){2-6} & 1 & 2 & 3 & 4 & 5 \\\\\n" \
	      TABLE_MID, f);
	for (u32 idx = 0; idx < LEN_LOGN; ++idx) {
		int logN = (int)params_logsN[idx];
		Params *p = params_single + pidx;
		f64 *tsgl = times + pidx;
		f64 *tfld = times + LEN_SINGLE + *idx0 + pidx;

		f64 x0 = tsgl[0] / tfld[0];
		f64 x1 = tsgl[1] / tfld[1];
		f64 x2 = tsgl[2] / tfld[2];
		f64 x3 = tsgl[3] / tfld[3];
		f64 x4 = tsgl[4] / tfld[4];
		#define Q6_FMT "    %2d & $%4.2f \\times$ & $%4.2f \\times$ & $%4.2f \\times$ & $%4.2f \\times$ & $%4.2f \\times$ \\\\\n"
		fprintf(f, Q6_FMT, logN, x0, x1, x2, x3, x4);

		u32 w = 0;
		while (p[w].N == 1U << logN)
			++w;
		pidx += w;
	}
	fputs(TABLE_END, f);

	fputs(TLONG(folding2, crrrrr)
	      "\\multicolumn{3}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "\\cmidrule(lr{1pt}){1-3}\\cmidrule(lr{1pt}){4-6} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & $\\omega$ & \\texttt{naïve} & \\texttt{folded} & Speedup \\\\\n"
	      "\\midrule\\endfirsthead\n"
	      "\\multicolumn{3}{c}{Parameters} & \\multicolumn{3}{c}{Time (ms)} \\\\\n"
	      "\\cmidrule(lr{1pt}){1-3}\\cmidrule(lr{1pt}){4-6} $\\log_2 N$ & \\multicolumn{1}{c}{$\\ell$} & $\\omega$ & \\texttt{naïve} & \\texttt{folded} & Speedup \\\\\n"
	      "\\midrule\\endhead\n", f);
	for (u32 idx = 0; idx < LEN_SINGLE; ++idx) {
		Params *p = params_single + idx;
		f64 *t = times + LEN_SINGLE + *idx0;
		int logN = (int)lround(log2((f64)p[0].N));
		f64 t0 = times[idx] * 1000.0, t1 = t[idx] * 1000.0;
		#define Q6_FMT_LONG "%2d & %2d & %2d & %6.1f & %6.1f & %4.2f \\\\\n"
		fprintf(f, Q6_FMT_LONG, logN, (int)p[0].L, (int)p[0].W, t0, t1, t0 / t1);
	}
	fputs(TLONG_END, f);
	*idx0 += LEN_SINGLE;
}

int
main(int argc, char **argv)
{
	f64 times[LEN_PARAMS] = { 0 };
	if (argc != 2)
		fprintf(stderr, "usage: %s <times.dat>\n", argv[0]), exit(1);
	FILE *fdat = fopen(argv[1], "r");
	if (fdat == 0)
		perror("fopen"), exit(1);
	assert(fread(times, sizeof *times, LEN_PARAMS, fdat) == LEN_PARAMS);
	fclose(fdat);

	u32 idx = 0;
	printq1(stdout, times, &idx);
	printq2(stdout, times, &idx);
	printq3(stdout, times, &idx);
	printq4(stdout, times, &idx);
	printq5(stdout, times, &idx);
	printq6(stdout, times, &idx);

	return 0;
}
