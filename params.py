#!/usr/bin/env python3
from math import ceil, floor, sqrt
from sympy import isprime


C="""
typedef struct {{ u64 *mods, *p, N, L, K, R, W, W2, mu, kappa; }} Params;
u64 params_p[1] = {{ 0x700001 }};
u64 params_logsN[] = {{ {} }};
u64 params_logsQP[] = {{ {} }};

{}

Params params_single[] = {{
    {}
}};
Params params_folded[] = {{
    {}
}};
Params params_double[] = {{
    {}
}};
Params params_switch[] = {{
    {}
}};
"""

LOGS = [ (14, 430), (15, 868), (16, 1747), (17, 3523) ]
FMT = "{{ params_{} + {:2d}, params_p, 1U << {}, {:2d}, {:2d}, {:2d}, {:2d}, {:2d}, 0, 0 }}"
FSW = "{{ params_{} + {:2d}, params_p, 1U << {}, {:2d}, {:2d}, {:2d}, {:2d}, {:2d}, {}, {} }}"


def genprimes(p, logN, num=1):
    diff = 1 << (logN + 1)
    primes = []
    while num > 0:
        p -= diff
        if isprime(p):
            primes.append(p)
            num -= 1
    return primes


def fmtprimes(name, primes):
    return f"u64 params_{name}[] = {{ {", ".join(map(lambda p: f"0x{p:016x}", primes))} }};"


with open("params.c", "w+") as f:
    iceil = lambda x: int(ceil(x))
    optW = lambda l, k: iceil(l / k)
    optW2 = lambda l: int(round(sqrt((l + 1) / 2 * l)))
    mods, single, folded, double, switch = [], [], [], [], []

    B60 = []
    for logN, logQP in LOGS:
        bits = 60
        LK = logQP // 60
        primes60 = genprimes((1 << 60) + 1, logN, LK)
        initLK = (1 << bits) + 1 if bits != 60 else primes60[-1]
        mods.append(fmtprimes(f"mod{logN}", genprimes(initLK, logN, LK) + primes60))
        mods.append(fmtprimes(f"mup{logN}", genprimes(initLK, logN + 1, 2 * LK)))
        B60.append((logN, LK))

    # 4.1  Is ω = 1 or ω = 2 better if we can choose (single-decomposition)?
    for logN, LK in B60:
        L = LK // 2
        folded.append(FMT.format(f"mod{logN}", 0, logN, L, L, 0, 1, 0))
        folded.append(FMT.format(f"mod{logN}", 0, logN, L, iceil(L / 2), 0, 2, 0))

    # 4.2  Can increasing N actually be worth it (single-decomposition)?
    for logN, LK in B60:
        L1, L2 = LK - 1, LK - 2
        folded.append(FMT.format(f"mod{logN}", 0, logN, L2, 2, 0, iceil(L2 / 2), 0))
        folded.append(FMT.format(f"mup{logN}", 0, logN + 1, L2, L2, 0, 1, 0))
        folded.append(FMT.format(f"mod{logN}", 0, logN, L1, 1, 0, L1, 0))
        folded.append(FMT.format(f"mup{logN}", 0, logN + 1, L1, L1, 0, 1, 0))

    # 4.3  Is the single- or the double-decomposition technique better?
    for logN, LK in B60:
        L = LK - 1
        W2 = optW2(L)
        R = iceil(LK / W2)
        for l in range(L, 0, -1):
            folded.append(FMT.format(f"mod{logN}", L - l, logN, l, 1, 0, l, 0))
            double.append(FMT.format(f"mod{logN}", L - l, logN, l, 1, R, l, W2))
        for l in range(L, 0, -1):
            w = optW(l, LK - l)
            k = iceil(l / w)
            folded.append(FMT.format(f"mod{logN}", LK - l - k, logN, l, k, 0, w, 0))
            w, k = l, 1
            w2 = optW2(l)
            r = iceil((l + 1) / w2)
            double.append(FMT.format(f"mod{logN}", LK - l - k, logN, l, k, r, w, w2))

    # 4.4  How large is the speedup using mostly large primes?
    # 4.5  How costly is replacing large with small primes?
    for logN, logQP in LOGS:
        LK = logQP // 54
        LK -= 1 - (LK % 2)
        L54 = LK - 1 # even
        L36 = L54 // 2 * 3
        primes = genprimes((1 << 36) + 1, logN, L36)
        primes += genprimes((1 << 54) + 1, logN, LK)
        mods.append(fmtprimes(f"msw{logN}", primes))
        for l54 in range(L54, 0, -2):
            l36 = l54 // 2 * 3
            w = optW(l54, LK - l54)
            k = iceil(l54 / w)
            folded.append(FMT.format(f"msw{logN}", L36 - l36, logN, l36, k, 0, w, 0))
            folded.append(FMT.format(f"msw{logN}", L36, logN, l54, k, 0, w, 0))
            switch.append(FSW.format(f"msw{logN}", L36 - 3, logN, 3 + l54, k, 0, w, 0, 2, 3))

    # 4.6  How large are the speedups from constant folding?
    for logN, LK in B60:
        for w in range(1, 11):
            k = LK // (w + 1)
            if k == 0:
                continue
            single.append(FMT.format(f"mod{logN}", 0, logN, w * k, k, 0, w, 0))
            folded.append(FMT.format(f"mod{logN}", 0, logN, w * k, k, 0, w, 0))

    f.write(C.format(
        ", ".join(map(lambda x: str(x[0]), LOGS)),
        ", ".join(map(lambda x: str(x[1]), LOGS)),
        "\n".join(mods),
        ",\n    ".join(single),
        ",\n    ".join(folded),
        ",\n    ".join(double),
        ",\n    ".join(switch),
    ))
