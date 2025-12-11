/* bedavg.c
 *
 * Minimal single-threaded C version of the provided C++ program.
 * Compile:
 *   gcc -std=c99 -O2 -o bedavg bedavg.c
 *
 * Usage:
 *   ./bedavg met.tsv regions.bed out.tsv [--bedzero|--1based] [--inclusive]
 *
 * Notes:
 * - met.tsv: tab-separated header: chr\tpos\tSample1\tSample2...
 * - regions.bed: chr\tstart\tend\t[name]
 * - By default bed is treated as 0-based half-open (start <= pos-1 < end).
 * - --1based will treat regions as 1-based inclusive (start <= pos <= end).
 * - --inclusive makes matching use start <= pos <= end regardless of --bedzero/--1based.
 *
 * This is intentionally single-threaded and uses simple dynamic arrays and qsort.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef long long ll;

typedef struct {
    char *chr;
    ll start;
    ll end;
    char *name;
    int index; /* position after sorting by chr,start */
    /* accumulators filled at runtime */
    double *sums;
    ll *counts;
} Interval;

/* --- Helpers for dynamic arrays --- */
typedef struct {
    Interval *a;
    size_t n, cap;
} IntervalArray;

static void ia_init(IntervalArray *ia) { ia->a = NULL; ia->n = ia->cap = 0; }
static void ia_push(IntervalArray *ia, Interval iv) {
    if (ia->n == ia->cap) {
        ia->cap = ia->cap ? ia->cap * 2 : 16;
        ia->a = realloc(ia->a, ia->cap * sizeof(Interval));
        if (!ia->a) { perror("realloc"); exit(1); }
    }
    ia->a[ia->n++] = iv;
}

/* string trim (leading/trailing whitespace) */
static char *strtrim(char *s) {
    if (!s) return NULL;
    while (isspace((unsigned char)*s)) s++;
    if (*s == 0) return strdup("");
    char *e = s + strlen(s) - 1;
    while (e > s && isspace((unsigned char)*e)) { *e = '\0'; e--; }
    return strdup(s);
}

/* fast split by tab: returns array of pointers (allocated) and count in *cnt; caller must free array (but not the strings) */
static char **split_tab(char *line, int *cnt) {
    int cap = 8;
    char **out = malloc(sizeof(char*) * cap);
    int n = 0;
    char *p = line;
    while (1) {
        char *q = p;
        while (*q && *q != '\t' && *q != '\r' && *q != '\n') q++;
        size_t len = q - p;
        char *cell = malloc(len + 1);
        memcpy(cell, p, len);
        cell[len] = '\0';
        if (n == cap) { cap *= 2; out = realloc(out, sizeof(char*) * cap); }
        out[n++] = cell;
        if (!*q || *q == '\r' || *q == '\n') break;
        p = q + 1;
    }
    *cnt = n;
    return out;
}

/* free array returned by split_tab */
static void free_split_tab(char **arr, int n) {
    for (int i = 0; i < n; ++i) free(arr[i]);
    free(arr);
}

/* qsort comparator: by chr then start then end */
static int cmp_interval(const void *A, const void *B) {
    const Interval *a = (const Interval*)A;
    const Interval *b = (const Interval*)B;
    int c = strcmp(a->chr, b->chr);
    if (c) return c;
    if (a->start < b->start) return -1;
    if (a->start > b->start) return 1;
    if (a->end < b->end) return -1;
    if (a->end > b->end) return 1;
    return 0;
}

/* load bed file into IntervalArray */
static void read_bed(const char *bedfile, IntervalArray *ia) {
    FILE *f = fopen(bedfile, "r");
    if (!f) { perror(bedfile); exit(1); }
    char *line = NULL;
    size_t llen = 0;
    ssize_t l;
    while ((l = getline(&line, &llen, f)) != -1) {
        if (l == 0) continue;
        if (line[0] == '#') continue;
        /* split by tab safely */
        int nfields = 0;
        char *tmp = strdup(line);
        char **fields = split_tab(tmp, &nfields);
        if (nfields >= 3) {
            char *chr = strtrim(fields[0]);
            char *sstart = strtrim(fields[1]);
            char *send = strtrim(fields[2]);
            ll start = atoll(sstart);
            ll end = atoll(send);
            char *name = NULL;
            if (nfields >= 4) name = strtrim(fields[3]);
            else {
                /* generate name chr-start-end */
                size_t bufcap = strlen(chr) + 1 + 32 + 1 + 32 + 1;
                name = malloc(bufcap);
                snprintf(name, bufcap, "%s-%lld-%lld", chr, (long long)start, (long long)end);
            }
            Interval iv;
            iv.chr = chr;
            iv.start = start;
            iv.end = end;
            iv.name = name;
            iv.sums = NULL;
            iv.counts = NULL;
            iv.index = -1;
            ia_push(ia, iv);
            free(sstart); free(send);
        }
        free_split_tab(fields, nfields);
        free(tmp);
    }
    free(line);
    fclose(f);
}

/* parse met header line to get sample names; returns sample count and an array of strdup'd sample names (caller free) */
static char **parse_met_header(const char *header_line, int *n_samples) {
    char *tmp = strdup(header_line);
    int nfields = 0;
    char **fields = split_tab(tmp, &nfields);
    if (nfields < 3) {
        fprintf(stderr, "met header must have at least chr,pos,sample...\n");
        exit(1);
    }
    int ns = nfields - 2;
    char **samples = malloc(sizeof(char*) * ns);
    for (int i = 0; i < ns; ++i) samples[i] = strtrim(fields[2 + i]);
    *n_samples = ns;
    free_split_tab(fields, nfields);
    free(tmp);
    return samples;
}

/* main processing: walk through met file (assumed sorted). Use active interval approach.
 * We assume intervals array ia->a is sorted by chr,start.
 */
int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s met.tsv regions.bed out.tsv [--bedzero|--1based] [--inclusive]\n", argv[0]);
        return 1;
    }
    const char *metfile = argv[1];
    const char *bedfile = argv[2];
    const char *outfile = argv[3];
    int bed_zero_based = 1;
    int inclusive = 0;
    for (int i = 4; i < argc; ++i) {
        if (strcmp(argv[i], "--1based") == 0) bed_zero_based = 0;
        else if (strcmp(argv[i], "--bedzero") == 0) bed_zero_based = 1;
        else if (strcmp(argv[i], "--inclusive") == 0) inclusive = 1;
        else {
            /* ignore numeric threads param in minimal C version */
        }
    }

    IntervalArray ia;
    ia_init(&ia);
    read_bed(bedfile, &ia);
    if (ia.n == 0) { fprintf(stderr, "No intervals read from %s\n", bedfile); return 1; }

    /* sort intervals by chr and start */
    qsort(ia.a, ia.n, sizeof(Interval), cmp_interval);
    for (size_t i = 0; i < ia.n; ++i) ia.a[i].index = (int)i;

    /* open met file and read header */
    FILE *mf = fopen(metfile, "r");
    if (!mf) { perror(metfile); return 1; }
    char *mline = NULL;
    size_t mlen = 0;
    ssize_t mread;

    if ((mread = getline(&mline, &mlen, mf)) == -1) {
        fprintf(stderr, "Empty met file or can't read header\n");
        return 1;
    }
    /* remove trailing newline */
    while (mread > 0 && (mline[mread-1] == '\n' || mline[mread-1] == '\r')) { mline[mread-1] = '\0'; --mread; }

    int nsamples = 0;
    char **sample_names = parse_met_header(mline, &nsamples);

    /* allocate accumulators per interval */
    for (size_t i = 0; i < ia.n; ++i) {
        ia.a[i].sums = calloc((size_t)nsamples, sizeof(double));
        ia.a[i].counts = calloc((size_t)nsamples, sizeof(ll));
    }

    /* open output and write header */
    FILE *of = fopen(outfile, "w");
    if (!of) { perror(outfile); return 1; }
    fprintf(of, "name");
    for (int i = 0; i < nsamples; ++i) fprintf(of, "\t%s", sample_names[i]);
    fprintf(of, "\n");

    /* PROCESS: read met line by line; maintain active intervals for current chromosome */
    /* We'll iterate intervals using an index pointer into ia.a */
    size_t interval_ptr = 0; /* first interval not yet activated for current chr */
    size_t active_cap = 64;
    size_t *active = malloc(active_cap * sizeof(size_t));
    size_t active_n = 0;

    char *linebuf = NULL;
    size_t linecap = 0;
    char *current_chr = NULL;

    while ((mread = getline(&linebuf, &linecap, mf)) != -1) {
        if (mread <= 0) continue;
        /* trim newline */
        while (mread > 0 && (linebuf[mread-1] == '\n' || linebuf[mread-1] == '\r')) { linebuf[mread-1] = '\0'; --mread; }
        if (linebuf[0] == '\0') continue;

        /* split the met row quickly */
        int ncols = 0;
        char *row = strdup(linebuf);
        char **cols = split_tab(row, &ncols);
        if (ncols < 2) { free_split_tab(cols,ncols); free(row); continue; }
        char *rchr = strtrim(cols[0]);
        char *posstr = strtrim(cols[1]);
        ll pos = atoll(posstr);

        /* if chromosome changed, clear active and advance interval_ptr to first interval on this chr */
        if (!current_chr || strcmp(current_chr, rchr) != 0) {
            /* flush active */
            active_n = 0;
            /* find first interval with chr == rchr */
            size_t i = 0;
            while (i < ia.n && strcmp(ia.a[i].chr, rchr) < 0) ++i;
            /* now i is first >= rchr; if equal, set interval_ptr = i, else there are no intervals for this chr */
            if (i < ia.n && strcmp(ia.a[i].chr, rchr) == 0) interval_ptr = i;
            else interval_ptr = i; /* points to first interval after chr */
            /* replace current_chr */
            free(current_chr);
            current_chr = rchr; /* note: rchr is strdup'd by strtrim; we want to own it */
            /* free some temporaries */
            free(posstr);
            free_split_tab(cols, ncols);
            free(row);
            continue; /* we will re-process this line after resetting? No — we skipped adding active intervals; better to re-read: simpler: rewind file pointer back by line and continue outer loop doesn't support re-reading -> instead process this row now: we already advanced interval_ptr appropriate for this chr; so continue processing below using current_chr==rchr */
        }

        /* Activate intervals whose start is <= pos (adjust criterion based on bed_zero_based/inclusive) */
        while (interval_ptr < ia.n && strcmp(ia.a[interval_ptr].chr, rchr) == 0) {
            int should_add = 0;
            if (inclusive) {
                if (ia.a[interval_ptr].start <= pos) should_add = 1;
            } else if (bed_zero_based) {
                /* interval applies if start <= pos-1 (conservative) */
                if (ia.a[interval_ptr].start <= pos - 1) should_add = 1;
            } else {
                if (ia.a[interval_ptr].start <= pos) should_add = 1;
            }
            if (should_add) {
                if (active_n == active_cap) {
                    active_cap *= 2;
                    active = realloc(active, active_cap * sizeof(size_t));
                    if (!active) { perror("realloc active"); exit(1); }
                }
                active[active_n++] = interval_ptr;
                interval_ptr++;
            } else break;
        }

        /* remove expired from active */
        size_t out_i = 0;
        for (size_t ai = 0; ai < active_n; ++ai) {
            size_t idx = active[ai];
            ll istart = ia.a[idx].start;
            ll iend = ia.a[idx].end;
            int still = 0;
            if (inclusive) {
                if (istart <= pos && pos <= iend) still = 1;
            } else if (bed_zero_based) {
                ll p0 = pos - 1;
                if (istart <= p0 && p0 < iend) still = 1;
            } else {
                if (istart <= pos && pos <= iend) still = 1;
            }
            if (still) active[out_i++] = idx;
        }
        active_n = out_i;

        if (active_n == 0) { /* no interval covering this pos */ free(posstr); free_split_tab(cols,ncols); free(row); continue; }

        /* check that we have enough columns for samples */
        if (ncols < 2 + nsamples) {
            free(posstr); free_split_tab(cols,ncols); free(row); continue;
        }

        /* for each active interval, add sample values */
        for (size_t k = 0; k < active_n; ++k) {
            size_t idx = active[k];
            Interval *iv = &ia.a[idx];
            for (int s = 0; s < nsamples; ++s) {
                char *sv = cols[2 + s];
                if (!sv || sv[0] == '\0' || strcmp(sv, ".") == 0 || strcmp(sv, "-") == 0 || strcmp(sv, "NA") == 0) continue;
                /* use strtod for parsing */
                char *endptr = NULL;
                double v = strtod(sv, &endptr);
                if (endptr == sv) continue; /* not a number */
                iv->sums[s] += v;
                iv->counts[s] += 1;
            }
        }

        free(posstr);
        free_split_tab(cols, ncols);
        free(row);
    } /* end reading met */

    /* write results */
    for (size_t i = 0; i < ia.n; ++i) {
        fprintf(of, "%s", ia.a[i].name);
        for (int s = 0; s < nsamples; ++s) {
            fprintf(of, "\t");
            if (ia.a[i].counts == 0) fprintf(of, ".");
            else {
                double avg = ia.a[i].sums[s] / (double)ia.a[i].counts[s];
                fprintf(of, "%.6f", avg);
            }
        }
        fprintf(of, "\n");
    }

    /* cleanup */
    fclose(mf);
    fclose(of);
    for (int i = 0; i < nsamples; ++i) free(sample_names[i]);
    free(sample_names);
    free(active);
    free(mline);
    free(linebuf);
    free(current_chr); /* may be NULL */

    for (size_t i = 0; i < ia.n; ++i) {
        free(ia.a[i].chr);
        free(ia.a[i].name);
        free(ia.a[i].sums);
        free(ia.a[i].counts);
    }
    free(ia.a);

    return 0;
}
