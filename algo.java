/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package newpackage;

import java.util.Arrays;

/**
 *
 * @author Administrator
 */
public class algo {

    //brute force
    public static void BF(char[] x, char[] y) {
        int i, j, m = x.length, n = y.length;

        for (j = 0; j <= n - m; ++j) {
            if (cmp(x, Arrays.copyOfRange(y, j, j+m))) {
                System.out.println(j);
            }
        }
    }

    // knuth-morris-partt
    public static int[] preKmp(char[] x) {
        int[] kmpNext = new int[x.length + 1];
        int i = 0;
        int j = kmpNext[0] = -1;

        while (i < x.length) {
            while (j > -1 && x[i] != x[j]) {
                j = kmpNext[j];
            }
            i++;
            j++;
            if (i == x.length - 1 || j == x.length - 1) {
                break;
            }
            if (x[i] == x[j]) {
                kmpNext[i] = kmpNext[j];
            } else {
                kmpNext[i] = j;
            }
        }
        return kmpNext;
    }

    public static void KMP(char[] x, char[] y) {
        int[] kmpNext = preKmp(x);
        int i = 0;
        int j = 0;
        while (j < y.length) {
            while (i > -1 && x[i] != y[j]) {
                i = kmpNext[i];
            }
            i++;
            j++;
            if (i >= x.length - 1) {
                System.out.println(j - i);
                i = kmpNext[i];
            }
        }
    }

    // morris-partt
    public static int[] preMP(char[] x) {
        int[] mpNext = new int[x.length + 1];
        int i = 0;
        int j = mpNext[0] = -1;
        while (i < x.length) {
            while (j > -1 && x[i] != x[j]) {
                j = mpNext[j];
            }
            mpNext[++i] = ++j;
        }
        return mpNext;
    }

    public static void MP(char[] x, char[] y) {
        int[] mpNext = preMP(x);
        int i = 0;
        int j = 0;
        while (j < y.length) {
            while (i > -1 && x[i] != y[j]) {
                i = mpNext[i];
            }
            i++;
            j++;
            if (i >= x.length - 1) {
                System.out.println(j - i);
                i = mpNext[i];
            }
        }
    }

    // karp-rabin
    public static int hash(char[] x) {
        int result = 0;
        for (int i = 0; i < x.length; i++) {
            result = (int) (result + (int) x[i] * Math.pow(2, x.length - i - 1));
        }
        return result;
    }
    // tính giá trị băm của 1 mảng kí tự dựa t rên giá trị băm đã tính trước đó
    // a: là kí tự đầu tiên của mảng trước đó
    // b: là kí tự cuối cùng của mảng mới
    // h: là giá trị băm của mảng cũ
    // m: là chiều dài của chuỗi cần so sánh

    public static int rehash(char a, char b, int h, int m) {
        int result = 0;
        result = (int) ((2 * (h - a * Math.pow(2, m - 1))) + b);
        return result;
    }
    // in ra các vị trí của x trong y

    public static void KR(char[] x, char[] y) {
        int hx = hash(x);
        char[] init = Arrays.copyOfRange(y, 0, x.length);
        int hy = hash(init);
        if (hx == hy) {
            System.out.println("0");
        }
        for (int i = 1; i <= y.length - x.length; i++) {
            hy = rehash(y[i - 1], y[i + x.length - 1], hy, x.length);
            if (hx == hy) {
                System.out.println(i);
            }
        }
    }

    // not so naive
    public static boolean cmp(char[] x, char[] y) {
        if (x.length != y.length) {
            return false;
        }
        for (int i = 0; i < x.length; i++) {
            if (x[i] != y[i]) {
                return false;
            }
        }
        return true;
    }

    public static void NSN(char[] x, char[] y) {
        int k = 1;
        int l = 2;
        char[] x1 = Arrays.copyOfRange(x, 2, x.length);
        if (x[0] == x[1]) {
            k = 2;
            l = 1;
        }
        int i = 0;
        while (i <= y.length - x.length) {
            if (x[1] != y[i + 1]) {
                i = i + k;
            } else {
                if (x[0] == y[i] && cmp(x1, Arrays.copyOfRange(y, i + 2, i + x.length))) {
                    System.out.println(i);
                }
                i += l;
            }
        }
    }

    // boyer moore
    public static int[] preBmBc(char[] x) {
        int i, m = x.length;
        int[] bmBc = new int[256];
        for (i = 0; i < 256; i++) {
            bmBc[i] = m;
        }
        for (i = 0; i < m - 1; i++) {
            bmBc[x[i]] = m - i - 1;
        }
        return bmBc;
    }

    public static int[] suffixes(char[] x) {
        int i, f, g, m = x.length;
        int[] suff = new int[256];
        suff[m - 1] = m;
        f = m;
        g = m - 1;
        for (i = m - 2; i >= 0; i--) {
            if (i > g && suff[i + m - 1 - f] < i - g) {
                suff[i] = suff[i + m - 1 - f];
            } else {
                if (i < g) {
                    g = i;
                }
                f = i;
                while (g >= 0 && x[g] == x[g + m - 1 - f]) {
                    --g;
                }
                suff[i] = f - g;
            }
        }
        return suff;
    }

    public static int[] preBmGs(char[] x) {
        int i, j, m = x.length;
        int[] bmGs = new int[256];
        int[] suff = suffixes(x);
        for (i = 0; i < m; i++) {
            bmGs[i] = m;
        }
        j = 0;
        for (i = m - 1; i >= -1; i--) {
            if (i == -1 || suff[i] == i + 1) {
                for (; j < m - 1 - i; j++) {
                    if (bmGs[j] == m) {
                        bmGs[j] = m - 1 - i;
                    }
                }
            }
        }
        for (i = 0; i <= m - 2; i++) {
            bmGs[m - 1 - suff[i]] = m - 1 - i;
        }
        return bmGs;
    }

    public static void BM(char[] x, char[] y) {
        int i, j, m = x.length, n = y.length;
        int[] bmGs = preBmGs(x);
        int[] bmBc = preBmBc(x);

        j = 0;
        while (j <= n - m) {
            i = m - 1;
            while (i >= 0 && x[i] == y[i + j]) {
                i--;
            }
            if (i < 0) {
                System.out.println(j);
                j += bmGs[0];
            } else {
                j += Math.max(bmGs[i], bmBc[y[i + j]] - m + 1 + i);
            }
        }
    }

    // Turbo BM
    public static void TBM(char[] x, char[] y) {
        int i, j, bcShift, shift, u, v, turboShift, m = x.length, n = y.length;
        int[] bmGs = preBmGs(x);
        int[] bmBc = preBmBc(x);

        j = u = 0;
        shift = m;
        while (j <= n - m) {
            i = m - 1;
            while (i >= 0 && x[i] == y[i + j]) {
                i--;
                if (u != 0 && i == m - 1 - shift) {
                    i -= u;
                }
            }
            if (i < 0) {
                System.out.println(j);
                shift = bmGs[0];
                u = m - shift;
            } else {
                v = m - 1 - i;
                turboShift = u - v;
                bcShift = bmBc[y[i + j]] - m + 1 + i;
                shift = Math.max(turboShift, bcShift);
                shift = Math.max(shift, bmGs[i]);
                if (shift == bmGs[i]) {
                    u = Math.min(m - shift, v);
                } else {
                    if (turboShift < bcShift) {
                        shift = Math.max(shift, u + 1);
                    }
                    u = 0;
                }
            }
            j += shift;
        }
    }

    // horspool
    public static void HORSPOOL(char[] x, char[] y) {
        int j, m = x.length, n = y.length;
        int[] bmBc = preBmBc(x);
        char c;

        j = 0;
        while (j <= n - m) {
            c = y[j + m - 1];
            if (x[m - 1] == c && cmp(x, Arrays.copyOfRange(y, j, j + m))) {
                System.out.println(j);
            }
            j += bmBc[c];
        }
    }

    // quick search
    public static int[] preQsBc(char[] x) {
        int i, m = x.length;
        int[] qsBc = new int[256];
        for (i = 0; i < 256; i++) {
            qsBc[i] = m + 1;
        }
        for (i = 0; i < m; i++) {
            qsBc[x[i]] = m - i;
        }
        return qsBc;
    }

    public static void QS(char[] x, char[] y) {
        int j, m = x.length, n = y.length;
        int[] qsBc = preQsBc(x);

        j = 0;
        while (j <= n - m) {
            if (cmp(x, Arrays.copyOfRange(y, j, j + m))) {
                System.out.println(j);
            }
            if (j != n - m) {
                j += qsBc[y[j + m]];
            } else {
                break;
            }
        }
    }

    // tuned BM
    public static void TUNEDBM(char[] x, char[] y) {
        int j, k, shift, m = x.length, n = y.length;
        int[] bmBc = preBmBc(x);
        shift = bmBc[x[m - 1]];
        bmBc[x[m - 1]] = 0;

        j = 0;
        while (j <= n - m) {
            k = bmBc[y[j + m - 1]];
            while (k != 0) {
                j += k;
                k = bmBc[y[j + m - 1]];
                j += k;
                k = bmBc[y[j + m - 1]];
                j += k;
                k = bmBc[y[j + m - 1]];
            }
            if (cmp(x, Arrays.copyOfRange(y, j, j + m))) {
                System.out.println(j);
            }
            j += shift;
        }
    }

    // ZT
    public static int[][] preZtBc(char[] x) {
        int[][] ztBc = new int[256][256];
        int i, j, m = x.length;
        for (i = 0; i < 256; i++) {
            for (j = 0; j < 256; j++) {
                ztBc[i][j] = m;
            }
        }
        for (i = 0; i < 256; i++) {
            ztBc[i][x[0]] = m - 1;
        }
        for (i = 1; i < m - 1; i++) {
            ztBc[x[i - 1]][x[i]] = m - 1 - i;
        }
        return ztBc;
    }

    public static void ZT(char[] x, char[] y) {
        int i, j, m = x.length, n = y.length;
        int[][] ztBc = preZtBc(x);
        int[] bmGs = preBmGs(x);

        j = 0;
        while (j <= n - m) {
            i = m - 1;
            while (i >= 0 && x[i] == y[i + j]) {
                i--;
            }
            if (i < 0) {
                System.out.println(j);
                j += bmGs[0];
            } else {
                j += Math.max(bmGs[i], ztBc[y[j + m - 2]][y[j + m - 1]]);
            }
        }
    }

    //BR
    public static int[][] preBrBc(char[] x) {
        int a, b, i, m = x.length;
        int[][] brBc = new int[256][256];

        for (a = 0; a < 256; a++) {
            for (b = 0; b < 256; b++) {
                brBc[a][b] = m + 2;
            }
        }
        for (a = 0; a < 256; a++) {
            brBc[a][x[0]] = m + 1;
        }
        for (i = 0; i < m - 1; i++) {
            brBc[x[i]][x[i + 1]] = m - i;
        }
        for (a = 0; a < 256; a++) {
            brBc[x[m - 1]][a] = 1;
        }

        return brBc;
    }

    public static void BR(char[] x, char[] y) {
        int j, m = x.length, n = y.length;
        int[][] brBc = preBrBc(x);

        j = 0;

        while (j <= n - m) {
            if (cmp(x, Arrays.copyOfRange(y, j, j + m))) {
                System.out.println(j);
            }
            if (j < n - m - 1) {
                j += brBc[y[j + m]][y[j + m + 1]];
            } else {
                break;
            }
        }
    }

    // Smith
    public static void SMITH(char[] x, char[] y) {
        int j, m = x.length, n = y.length;
        int[] bmBc = preBmBc(x);
        int[] qsBc = preQsBc(x);

        j = 0;
        while (j <= n - m) {
            if (cmp(x, Arrays.copyOfRange(y, j, j + m))) {
                System.out.println(j);
            }
            if (j < n - m) {
                j += Math.max(bmBc[y[j + m - 1]], qsBc[y[j + m]]);
            } else {
                break;
            }
        }
    }

    // Raita
    public static void RAITA(char[] x, char[] y) {
        int j, m = x.length, n = y.length;
        int[] bmBc = preBmBc(x);

        j = 0;
        while (j <= n - m) {
            char c = y[j + m - 1];
            char first = x[0];
            char middle = x[m / 2];
            char last = x[m - 1];
            if (first == y[j] && middle == y[j + m / 2] && last == y[j + m - 1]
                    && cmp(Arrays.copyOfRange(x, 1, m - 1), Arrays.copyOfRange(y, j + 1, j + m - 1))) {
                System.out.println(j);
            }
            j = j + bmBc[c];
        }
    }

    // KMP Skip
    public static int attempt(char[] x, char[] y, int start, int wall) {
        int k = wall - start;
        while (k < x.length && x[k] == y[k + start]) {
            k++;
        }
        return k;
    }

    public static void KMPSKIP(char[] x, char[] y) {
        int m = x.length, n = y.length;
        int i, j, k, kmpStart, per, start, wall;
        int[] list = new int[m + 1];
        int[] z = new int[256];
        int[] kmpNext = preKmp(x);
        int[] mpNext = preMP(x);
        Arrays.fill(z, -1);
        Arrays.fill(list, -1);
        z[x[0]] = 0;
        for (i = 1; i < m; i++) {
            list[i] = z[x[i]];
            z[x[i]] = 1;
        }

        wall = 0;
        per = m - kmpNext[m - 1];
        i = j = -1;
        do {
            j += m;
        } while (j < n && z[y[j]] < 0);
        if (j >= n) {
            return;
        }
        i = z[y[j]];
        start = j - i - 1;

        while (start <= n - m) {
            if (start > wall) {
                wall = start;
            }
            k = attempt(x, y, start, wall);
            wall = start + k;
            if (k == m) {
                System.out.println(start);
                i -= per;
            } else {
                i = list[i];
            }
            if (i < 0) {
                do {
                    j += m;
                } while (j < n && z[y[j]] < 0);
                if (j >= n) {
                    return;
                }
                i = z[y[j]];
            }
            kmpStart = start + k + kmpNext[k];
            k = kmpNext[k];
            start = j - i;
            while (start < kmpStart || kmpStart < start && start < wall) {
                if (start < kmpStart) {
                    i = list[i];
                    if (i < 0) {
                        do {
                            j += m;
                        } while (j < n && z[y[j]] < 0);
                        if (j >= n) {
                            return;
                        }
                        i = z[y[j]];
                    }
                    start = j - i;
                } else {
                    kmpStart += k - mpNext[k];
                    k = mpNext[k];
                }
            }
        }
    }

    // AXAMAC
    public static void AXAMAC(char[] x, char[] y) {
        int i, j, k, ell, m = x.length, n = y.length;
        int[] kmpNext = preKmp(x);
        ell = 1;
        while (x[ell - 1] == x[ell]) {
            ell++;
        }
        if (ell == m) {
            ell = 0;
        }

        i = ell;
        j = k = 0;
        while (j <= n - m) {
            while (i < m && x[i] == y[i + j]) {
                i++;
            }
            if (i >= m) {
                while (k < ell && x[k] == y[j + k]) {
                    k++;
                }
                if (k >= ell) {
                    System.out.println(j);
                }
            }
            j += i - kmpNext[i];
            if (i == ell) {
                k = Math.max(0, k - 1);
            } else if (kmpNext[i] <= ell) {
                k = Math.max(0, kmpNext[i]);
                i = ell;
            } else {
                k = ell;
                i = kmpNext[i];
            }
        }
    }

    // colussi
    public static int preColussi(char[] x, int[] h, int[] shift, int next[]) {
        int i, k, nd, q, r = 0, s, m = x.length;
        int[] hmax = new int[m + 1];
        int[] kmin = new int[m + 1];
        int[] nhd0 = new int[m + 1];
        int[] rmin = new int[m + 1];

        /* Computation of hmax */
        i = k = 1;
        do {
            while (x[i] == x[i - k]) {
                i++;
            }
            hmax[k] = i;
            q = k + 1;
            while (hmax[q - k] + k < i) {
                hmax[q] = hmax[q - k] + k;
                q++;
            }
            k = q;
            if (k == i + 1) {
                i = k;
            }
            if (i == m) {
                break;
            }
        } while (k <= m);

        /* Computation of kmin */
        Arrays.fill(kmin, 0);
        for (i = m; i >= 1; --i) {
            if (hmax[i] < m) {
                kmin[hmax[i]] = i;
            }
        }

        /* Computation of rmin */
        for (i = m - 1; i >= 0; --i) {
            if (hmax[i + 1] == m) {
                r = i + 1;
            }
            if (kmin[i] == 0) {
                rmin[i] = r;
            } else {
                rmin[i] = 0;
            }
        }

        /* Computation of h */
        s = -1;
        r = m;
        for (i = 0; i < m; ++i) {
            if (kmin[i] == 0) {
                h[--r] = i;
            } else {
                h[++s] = i;
            }
        }
        nd = s;

        /* Computation of shift */
        for (i = 0; i <= nd; ++i) {
            shift[i] = kmin[h[i]];
        }
        for (i = nd + 1; i < m; ++i) {
            shift[i] = rmin[h[i]];
        }
        shift[m] = rmin[0];

        /* Computation of nhd0 */
        s = 0;
        for (i = 0; i < m; ++i) {
            nhd0[i] = s;
            if (kmin[i] > 0) {
                ++s;
            }
        }

        /* Computation of next */
        for (i = 0; i <= nd; ++i) {
            next[i] = nhd0[h[i] - kmin[h[i]]];
        }
        for (i = nd + 1; i < m; ++i) {
            next[i] = nhd0[m - rmin[h[i]]];
        }
        next[m] = nhd0[m - rmin[h[m - 1]]];

        return (nd);
    }

    public static void COLUSSI(char[] x, char[] y) {
        int i, j, last, nd, m = x.length, n = y.length;
        int[] h = new int[m + 1];
        int[] next = new int[m + 1];
        int[] shift = new int[m + 1];

        /* Processing */
        nd = preColussi(x, h, shift, next);

        /* Searching */
        i = j = 0;
        last = -1;
        while (j <= n - m) {
            while (i < m && last < j + h[i]
                    && x[h[i]] == y[j + h[i]]) {
                i++;
            }
            if (i >= m || last >= j + h[i]) {
                System.out.println(j);
                i = m;
            }
            if (i > nd) {
                last = j + m - 1;
            }
            j += shift[i];
            i = next[i];
        }
    }

    public static void main(String[] args) {
        String x = "GCAGAGAG";
        char[] X = x.toCharArray();
        String y = "GCATCGCAGAGAGTATACAGTACG";
        char[] Y = y.toCharArray();
        System.out.println("BF");
        BF(X,Y);
        System.out.println("KMP");
        KMP(X,Y);
        System.out.println("MP");
        MP(X,Y);
        System.out.println("KR");
        KR(X,Y);
        System.out.println("NSN");
        NSN(X,Y);
        System.out.println("BM");
        BM(X,Y);
        System.out.println("TBM");
        TBM(X,Y);
        System.out.println("HORSPOOL");
        HORSPOOL(X,Y);
        System.out.println("QS");
        QS(X,Y);
        System.out.println("TUNEDBM");
        TUNEDBM(X,Y);
        System.out.println("ZT");
        ZT(X,Y);
        System.out.println("BR");
        BR(X,Y);
        System.out.println("SMITH");
        SMITH(X,Y);
        System.out.println("RAITA");
        RAITA(X,Y);
        System.out.println("KMPSKIP");
        KMPSKIP(X,Y);
        System.out.println("AXAMAC");
        AXAMAC(X,Y);
//        COLUSSI(X,Y);
    }
}
