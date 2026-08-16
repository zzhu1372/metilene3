#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct {
    int64_t start, end;
} Interval;

typedef struct {
    char *chr;
    Interval *iv;
    size_t n, cap;
} Chrom;

static Chrom *chroms = NULL;
static size_t nchrom = 0, chrom_cap = 0;

static Chrom *get_chrom(const char *name) {
    for (size_t i = 0; i < nchrom; i++) {
        if (strcmp(chroms[i].chr, name) == 0)
            return &chroms[i];
    }

    if (nchrom == chrom_cap) {
        chrom_cap = chrom_cap ? chrom_cap * 2 : 32;
        chroms = realloc(chroms, chrom_cap * sizeof(*chroms));
        if (!chroms) {
            perror("realloc");
            exit(EXIT_FAILURE);
        }
    }

    Chrom *c = &chroms[nchrom++];
    c->chr = strdup(name);
    c->iv = NULL;
    c->n = 0;
    c->cap = 0;

    return c;
}

static void add_interval(Chrom *c, int64_t start, int64_t end) {
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 1024;
        c->iv = realloc(c->iv, c->cap * sizeof(*c->iv));

        if (!c->iv) {
            perror("realloc");
            exit(EXIT_FAILURE);
        }
    }

    c->iv[c->n++] = (Interval){start, end};
}

static int cmp_iv(const void *a, const void *b) {
    const Interval *x = a;
    const Interval *y = b;

    if (x->start < y->start) return -1;
    if (x->start > y->start) return 1;
    if (x->end < y->end) return -1;
    if (x->end > y->end) return 1;

    return 0;
}

static void merge_intervals(Chrom *c) {
    if (c->n == 0)
        return;

    qsort(c->iv, c->n, sizeof(*c->iv), cmp_iv);

    size_t out = 0;

    for (size_t i = 1; i < c->n; i++) {
        if (c->iv[i].start <= c->iv[out].end) {
            if (c->iv[i].end > c->iv[out].end)
                c->iv[out].end = c->iv[i].end;
        } else {
            c->iv[++out] = c->iv[i];
        }
    }

    c->n = out + 1;
}

static Chrom *find_chrom(const char *name) {
    for (size_t i = 0; i < nchrom; i++) {
        if (strcmp(chroms[i].chr, name) == 0)
            return &chroms[i];
    }

    return NULL;
}

static int covered(const Chrom *c, int64_t p) {
    size_t lo = 0;
    size_t hi = c->n;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (c->iv[mid].start <= p)
            lo = mid + 1;
        else
            hi = mid;
    }

    if (lo == 0)
        return 0;

    return p < c->iv[lo - 1].end;
}

static void load_bed(const char *path) {
    FILE *f = fopen(path, "r");

    if (!f) {
        perror(path);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t cap = 0;

    while (getline(&line, &cap, f) != -1) {
        if (line[0] == '#' || line[0] == '\n')
            continue;

        char *save = NULL;

        char *chr = strtok_r(line, "\t ", &save);
        char *s   = strtok_r(NULL, "\t ", &save);
        char *e   = strtok_r(NULL, "\t ", &save);

        if (!chr || !s || !e)
            continue;

        int64_t start = strtoll(s, NULL, 10);
        int64_t end   = strtoll(e, NULL, 10);

        if (end > start)
            add_interval(get_chrom(chr), start, end);
    }

    free(line);
    fclose(f);

    for (size_t i = 0; i < nchrom; i++)
        merge_intervals(&chroms[i]);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr,
                "Usage: %s met.tsv regions.bed out.tsv\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    const char *met_path = argv[1];
    const char *bed_path = argv[2];
    const char *out_path = argv[3];

    load_bed(bed_path);

    FILE *f = fopen(met_path, "r");

    if (!f) {
        perror(met_path);
        return EXIT_FAILURE;
    }

    FILE *out = fopen(out_path, "w");

    if (!out) {
        perror(out_path);
        fclose(f);
        return EXIT_FAILURE;
    }

    char *line = NULL;
    size_t cap = 0;

    /* Copy header */
    if (getline(&line, &cap, f) != -1)
        fputs(line, out);

    while (getline(&line, &cap, f) != -1) {

        char *tab1 = strchr(line, '\t');
        if (!tab1)
            continue;

        char *tab2 = strchr(tab1 + 1, '\t');
        if (!tab2)
            continue;

        /*
         * Temporarily terminate chr and pos.
         * We never parse the sample columns.
         */
        *tab1 = '\0';
        *tab2 = '\0';

        char *chr = line;
        int64_t pos = strtoll(tab1 + 1, NULL, 10);

        Chrom *c = find_chrom(chr);

        *tab1 = '\t';
        *tab2 = '\t';

        /*
         * Assumption:
         *   met.tsv pos = 1-based
         *   BED = 0-based half-open
         *
         * Therefore genomic position pos corresponds
         * to BED coordinate pos - 1.
         */
        if (c && covered(c, pos - 1))
            fputs(line, out);
    }

    free(line);

    fclose(f);
    fclose(out);

    for (size_t i = 0; i < nchrom; i++) {
        free(chroms[i].chr);
        free(chroms[i].iv);
    }

    free(chroms);

    return EXIT_SUCCESS;
}