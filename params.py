#!/usr/bin/env python3
from sage.all import *
import json

# commands:
#   for i in $(seq  0 33); do taskset -c 0 ./build/single${i} --benchmark_repetitions=${reps} --benchmark_out=data/single${i}.json && taskset -c 0 ./build/folded${i} --benchmark_repetitions=${reps} --benchmark_out=data/folded${i}.json; done && for i in $(seq 34 58); do taskset -c 0 ./build/double${i} --benchmark_repetitions=${reps} --benchmark_out=data/double${i}.json; done && for i in $(seq 59 82); do taskset -c 0 ./build/switch${i} --benchmark_repetitions=${reps} --benchmark_out=data/switch${i}.json; done
#   for i in $(seq  0 33); do scp ${server}:owl/data/single${i}.json data/single${i}.json; done
#   for i in $(seq  0 33); do scp ${server}:owl/data/folded${i}.json data/folded${i}.json; done
#   for i in $(seq 34 58); do scp ${server}:owl/data/double${i}.json data/double${i}.json; done
#   for i in $(seq 59 82); do scp ${server}:owl/data/switch${i}.json data/switch${i}.json; done
CMAKE = False
LATEX = True
REPS = 10

def bench2latex(params, single, folded, verbose=LATEX):
    if verbose:
        print((
            "    \\begin{tabular}{c*{5}{r}}\n"
            "    \\toprule\n"
            "    \\multicolumn{3}{c}{Parameters} & \\multicolumn{3}{c}{Time ($ms$)} \\\\\n"
            "    \\cmidrule(lr{1pt}){1-3}\\cmidrule(lr{1pt}){4-6} $\\log N$ & $\\log q$ & $\\omega$ & naïve & folded & speedup \\\\"))

    for idx, p, s, f in zip(range(len(params)), params, single, folded):
        logN, logq, omega = p[0], sum(p[3]), p[6]

        ts, tf = 0, 0
        for i in range(4):
            ts += s[i] / 4.0
            tf += f[i] / 4.0

        if verbose:
            if idx in [0, 6, 14, 24]:
                print("    \\midrule")
            print(f"    {logN} & {logq:4d} & {omega:2d} & {ts:6.1f} & {tf:6.1f} & {(ts / tf - 1) * 100:4.1f}\\,\\%\\\\")

    if verbose:
        print((
            "    \\bottomrule\n"
            "    \\end{tabular}\n"))


def double2latex(params, double, verbose=LATEX):
    if verbose:
        print((
            "    \\begin{tabular}{c*{4}{r}}\n"
            "    \\toprule\n"
            "    \\multicolumn{4}{c}{Parameters} & \\multicolumn{1}{c}{Time ($ms$)} \\\\\n"
            "    \\cmidrule(lr{1pt}){1-4} $\\log N$ & $\log q$ & $\\omega$ & $\\tilde\\omega$ \\\\"))

    for idx, p, d in zip(range(len(params)), params, double):
        logN, logq, omega, omega2 = p[0], sum(p[3]), p[6], p[8]
        t = 0
        for i in range(4):
            t += d[i] / 4.0

        if verbose:
            if idx in [0, 5, 11, 18]:
                print("    \\midrule")
            print(f"    {logN} & {logq:4d} & {omega:2d} & {omega2:2d} & {t:6.1f} \\\\")

    if verbose:
        print((
            "    \\bottomrule\n"
            "    \\end{tabular}\n"))


def est(logN, logq):
    N = 1 << logN
    alpha = 0.05
    beta = 0.33
    gamma = 17.88
    delta = 0.65

    return -log(alpha * logq / N, 2) * beta * N / logq + gamma * pow(logq / N, delta) * log(N / logq, 2)


def estdouble(l):
    est = []
    for w in range(1, l + 1):
        for w2 in range(1, l + 1):
            ntt = (w + 2 * w2) * (w * l + w2 * l + l) / (w * w2)
            est.append((ntt, w, w2))
    est.sort()

    _, w, w2 = est[0]
    assert w == l and w2 == estw2(l)

    est = []
    for w in range(1, l + 1):
        for w2 in range(1, l + 1):
            mul = l * (2 * w + 2 * w2 + 3 * l / w + l / w2 + l / (w * w2) + 5)
            est.append((mul, w, w2))
    est.sort()
    # TODO: closed formula?


def estmax(logN, verbose=False, eps=0):
    logq = 0
    while True:
        if est(logN, logq + 1) < 128 + eps:
            break
        logq += 1

    if verbose:
        print(f"[+] max{logN:02d}: {logq:5d}")
    return logq


def estw2(l):
    return round(sqrt(l * (l + 1) / 2))


def genprime(init, n, reverse=True):
    if reverse:
        p = init - (init % (2 * n)) - 2 * n + 1
        assert p < init and p % (2 * n) == 1
    else:
        p = init - (init % (2 * n)) + 2 * n + 1
        assert p > init and p % (2 * n) == 1

    while not is_prime(p):
        if reverse:
            p -= 2 * n
        else:
            p += 2 * n

    return p


def genparams(logN, b, beta=None, omega=2, omega2=1, logp=30, B=60, eps=0):
    N = 1 << logN
    logq = sum(b)
    logP = ceil(logq / omega)

    if not beta:
        k = ceil(logP / B)
        beta = [ceil(logP / k)] * k
    assert sum(beta) >= logP
    assert est(logN, sum(b) + sum(beta)) >= 128 - eps

    logP = sum(beta)
    logE = logN - 1 + logP + ceil((logq + logP) / omega2)
    r = ceil(logE / B)

    init = {}
    for bits in set(b + beta + [ceil(logE / r), logp]):
        init[bits] = 1 << bits

    p = genprime(init[logp], N)
    init[logp] = p

    q = []
    for bits in b:
        q.append(genprime(init[bits], N))
        init[bits] = q[-1]

    P = []
    for bits in beta:
        P.append(genprime(init[bits], N))
        init[bits] = P[-1]

    E = []
    bits = ceil(logE / r)
    for i in range(r):
        E.append(genprime(init[bits], N))
        init[bits] = E[-1]

    return logN, p, q, b, P, beta, omega, E, omega2


def params2cmake(idx, logN, p, q, b, P, beta, omega, E, omega2, switch=0, verbose=CMAKE):
    s = [
        str(1 << logN),
        f"\\\"{p}\\\"",
        ','.join(map(str, q)),
        ','.join(map(str, P)),
        str(omega),
        ','.join(map(str, E)),
        str(omega2),
        str(switch)
    ]

    if verbose:
        print(f"set(PARAM{idx} \"{';'.join(s)}\")")


def readbench(idx, params, keyswitch="", reps=REPS, offset=4, verbose=False):
    with open(f"data/{keyswitch}{idx}.json", "r") as f:
        data = json.load(f)
    mean = reps

    i = idx
    idx = mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_intt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    real0 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = reps + offset + mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_ntt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    real1 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 2 * (reps + offset) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_intt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    real2 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 3 * (reps + offset) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_ntt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    real3 = data["benchmarks"][idx]["real_time"] / 1000000

    n, b, k, omega, r, omega2 = params[0], params[3], len(params[4]), params[6], len(params[7]), params[8]
    if verbose:
        print(((
            f"idx: {i}, logN: {n}, logq: {sum(b)}, omega: {omega}\n"
            "==============================================\n"
            f"    domain        {keyswitch}\n"
            "---------------   --------\n"
            "input    output     (ms)\n"
            "--------------------------\n"
            f" coef     coef   {real0:8.2f}\n"
            f" coef     NTT    {real1:8.2f}\n"
            f" NTT      coef   {real2:8.2f}\n"
            f" NTT      NTT    {real3:8.2f}\n"
            "==============================================\n")))

    return [real0, real1, real2, real3]


def readswitch(idx, params, reps=REPS, verbose=False):
    with open(f"data/switch{idx}.json", "r") as f:
        data = json.load(f)
    mean = reps

    i = idx
    idx = mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_intt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    folded0 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = reps + 4 + mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_ntt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    folded1 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 2 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_intt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    folded2 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 3 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_ntt_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    folded3 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 4 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_intt_replace_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    replace0 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 5 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_intt_to_ntt_replace_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    replace1 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 6 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_intt_replace_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    replace2 = data["benchmarks"][idx]["real_time"] / 1000000

    idx = 7 * (reps + 4) + mean
    assert data["benchmarks"][idx]["name"] == f"BM_ntt_to_ntt_replace_mean"
    assert data["benchmarks"][idx]["time_unit"] == "ns"
    replace3 = data["benchmarks"][idx]["real_time"] / 1000000

    n, b, k, omega, r, omega2 = params[0], params[3], len(params[4]), params[6], len(params[7]), params[8]
    if verbose:
        print(((
            f"idx: {i}, logN: {n}, logq: {sum(b)}, omega: {omega}\n"
            "=====================================\n"
            f"    domain            switch (ms)\n"
            "---------------   -------------------\n"
            "input    output    folded    replace \n"
            "-------------------------------------\n"
            f" coef     coef    {folded0:8.2f}   {replace0:8.2f}\n"
            f" coef     NTT     {folded1:8.2f}   {replace1:8.2f}\n"
            f" NTT      coef    {folded2:8.2f}   {replace2:8.2f}\n"
            f" NTT      NTT     {folded3:8.2f}   {replace3:8.2f}\n"
            "=====================================\n")))

    return [(folded0, replace0), (folded1, replace1), (folded2, replace2), (folded3, replace3)]


def switch2latex(params, switch, verbose=LATEX):
    if verbose:
        print((
            "    \\begin{tabular}{c*{6}{r}}\n"
            "    \\toprule\n"
            "    $\\log N$ & $\\log q$ & \\multicolumn{2}{c}{$\\log q_i$} & \\multicolumn{3}{c}{Time ($ms$)} \\\\\n"
            "    \\cmidrule(lr{1pt}){3-4}\\cmidrule(lr{1pt}){5-7} && 36 & 54 & folded & switch & overhead \\\\"))

    for idx, p, s in zip(range(2, len(params)), params[2:], switch[2:]):
        logN, b, omega = p[0], p[3], p[6]

        tf, ts = 0, 0
        for i in range(4):
            tf += s[i][0] / 4.0
            ts += s[i][1] / 4.0

        if verbose:
            if idx in [2, 6, 12, 18]:
                print("    \\midrule")
            if b.count(54) > 0:
                print(f"                        &                       & {b.count(36):2d} & {b.count(54):2d} & {tf:6.1f} & {ts:6.1f} & {(ts / tf - 1) * 100:4.1f}\\,\\%\\\\")
            else:
                print(f"    \\multirow{{2}}{{*}}{{{logN}}} & \\multirow{{2}}{{*}}{{{sum(b):4d}}} & {b.count(36):2d} &    & {tf:6.1f} &        &         \\\\")

    if verbose:
        print((
            "    \\bottomrule\n"
            "    \\end{tabular}"))


if __name__ == "__main__":
    logN = [14, 15, 16, 17]
    B = 60

    for l in range(1, 201):
        estdouble(l)

    # print("[+] params, pdouble")
    params, pdouble = [], []
    for n in logN:
        logmax = estmax(n, eps=2)

        for factor in [0.5, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95]:
            logq = factor * logmax
            l = ceil(logq / B)
            b = [ceil(logq / l)] * l

            wmin = ceil(logq / (logmax - logq))
            if wmin > l:
                continue
            w2 = estw2(l)

            params.append(genparams(n, b, omega=wmin, B=B, eps=.5))
            pdouble.append(genparams(n, b, omega=l, omega2=w2, B=B, eps=.5))
            if factor == 0.5:
                params.append(genparams(n, b, omega=2, B=B, eps=.5))
                assert wmin == 1
            elif factor >= 0.9:
                params.append(genparams(n + 1, b, omega=2, B=B, eps=.5))

    # print("[+] pswitch")
    pswitch = []
    for n in logN:
        logmax = estmax(n, eps=2)

        for factor in [1/2, 2/3, 3/4]:
            l36 = floor(factor * logmax / 36)
            logq = l36 * 36
            assert l36 >= 4

            l54 = (l36 - 3) // 3 * 2
            b = [54] * l54 + [36] * (l36 % 3 + 3)
            assert sum(b) == logq

            w = ceil(logq / (logmax - logq))
            pswitch.append(genparams(n, [36] * l36, omega=w, B=B, eps=.5))
            pswitch.append(genparams(n, b, omega=w, B=B, eps=.5))

    # print("[+] cmake")
    idx = 0
    for p in params:
        params2cmake(idx, *p)
        idx += 1
    # print(f"[idx] {idx:2d}")
    assert idx == 34

    for p in pdouble:
        params2cmake(idx, *p)
        idx += 1
    # print(f"[idx] {idx:2d}")
    assert idx == 59

    for p in pswitch:
        params2cmake(idx, *p, switch=3)
        idx += 1
    # print(f"[idx] {idx:2d}")
    assert idx == 83

    idx = 0
    single, folded = [], []
    for p in params:
        single.append(readbench(idx, p, "single"))
        folded.append(readbench(idx, p, "folded"))
        idx += 1

    double = []
    for p in pdouble:
        double.append(readbench(idx, p, "double"))
        idx += 1

    switch = []
    for p in pswitch:
        switch.append(readswitch(idx, p))
        idx += 1

    # print("[+] latex")
    bench2latex(params, single, folded)
    double2latex(pdouble, double)
    switch2latex(pswitch, switch)
