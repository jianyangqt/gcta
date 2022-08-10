#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>

//#define ACAT_UTILITY

//When compile under windows, M_PI is defined in corecrt_math_defines.h
#if defined _WIN32 || defined _WIN64
#include <corecrt_math_defines.h>
#endif

#define BETA_DENSITY(x, a, b)  \
    (((gsl_sf_gamma(a + b)) / (gsl_sf_gamma(a) * gsl_sf_gamma(b))) * pow((x), a - 1.0) * pow(1.0 - (x), b - 1.0))

#define PCAUCHY(x) \
    ((1.0L / M_PI) * atan((x - 0L)/1.0L) + 0.5L)

#define PCAUCHY2(x) \
    atan(1.0 / (x)) / M_PI

#define HASH_LEN 100
#define LINE_BUFFER_LEN 4096
#define FIELD_BUFFER_LEN 512

#ifndef LOCAL_STRUCT /*LOCAL_STRUCT*/
#define LOCAL_STRUCT

struct SNP_NODE{
    unsigned long     snp_pos;
    double            weight;
    double            p_value;
    struct SNP_NODE * next;
};

struct GENE_NODE{
    unsigned long             gene_pos_s_e[2];
    char                      gene_name[20];
    struct SNP_POINTER_NODE * snp_pointer;
    struct SNP_POINTER_NODE * snp_pointer_tail;
    struct GENE_NODE        * prev;
    struct GENE_NODE        * next;
};

struct HASH_NODE{
    struct SNP_NODE  *   snp_dt;
    struct SNP_NODE  *   snp_dt_tail;
    struct GENE_NODE *   gene_dt;
    struct GENE_NODE *   gene_dt_tail;
    unsigned long    *   gene_order_list;
    unsigned long        gene_num;
};

struct SNP_POINTER_NODE{
    struct SNP_NODE         * snp_pt;
    struct SNP_POINTER_NODE * next;
};
#endif /*LOCAL_STRUCT*/

static unsigned int hash_func(char *);
static const char * unhash_func(unsigned int);
static void structure_gene_data(struct HASH_NODE *, const char *, bool, unsigned int);
static void structure_snp_data(struct HASH_NODE *, const char *, bool, double, unsigned int);
static bool judge_snp_in_gene(unsigned long, unsigned long *, unsigned long);
static void * thread_worker_1(void *);
static void print_res_file(struct HASH_NODE *, const char *);
static void free_mem(struct HASH_NODE *);

int acat_func(const char *, const char *, double, unsigned int, unsigned int, const char *);


#ifdef ACAT_UTILITY
int
main (int argc, char * argv[])
{
    const char * gene_file_name;
    const char * snp_file_name;
    double max_af = 0.01;
    unsigned int min_sample = 10;
    unsigned int expand_len = 0;
    const char * res_file_name = "test_result.csv";

    if (argc != 3){
        fprintf(stderr, "error, argument number error.\n");
        exit(EXIT_FAILURE);
    }

    gene_file_name = argv[1];
    snp_file_name = argv[2];

    acat_func(gene_file_name, snp_file_name, max_af, min_sample, expand_len, res_file_name);
    return 0;
}
#endif


int
acat_func (const char * gene_list_file, const char * snp_list_file, double max_af,\
    unsigned int min_sample, unsigned int extend_len, const char * out_f)
{

    printf("args: --acat --gene-list %s --snp-list %s --max-maf %f --min-mac  %u --wind %u --out %s\n", \
        gene_list_file, snp_list_file, max_af, min_sample, extend_len, out_f);
    unsigned char i = 0;
    time_t t_start, t_stop;

    struct HASH_NODE hash_dt[HASH_LEN];
    for (i = 0; i < HASH_LEN; ++i){
        hash_dt[i].snp_dt = NULL;
        hash_dt[i].snp_dt_tail = NULL;
        hash_dt[i].gene_dt = NULL;
        hash_dt[i].gene_dt_tail = NULL;
        hash_dt[i].gene_order_list = NULL;
        hash_dt[i].gene_num = 0;
    }


    printf(">structure gene_list data start...\n");
    time(&t_start);
    structure_gene_data(hash_dt, gene_list_file, false, extend_len);
    time(&t_stop);
    printf("<down. used %lf second(s)\n\n", difftime(t_stop, t_start));


    printf(">structure snp_list data start...\n");
    time(&t_start);
    structure_snp_data(hash_dt, snp_list_file, true, max_af, min_sample);
    time(&t_stop);
    printf("<down. used %lf second(s)\n\n", difftime(t_stop, t_start));


    printf(">mapping snp to gene start...\n");
    time(&t_start);
    for (i = 0; i < HASH_LEN; i++){
        if ( hash_dt[i].gene_dt != NULL && hash_dt[i].snp_dt != NULL ){
            thread_worker_1(&(hash_dt[i]));
        }
    }
    time(&t_stop);
    printf("<down. used %lf second(s)\n\n", difftime(t_stop, t_start));
    

    printf(">calculate cauchy and print results start...\n");
    time(&t_start);
    print_res_file(hash_dt, out_f);
    time(&t_stop);
    printf("<down. used %lf second(s)\n\n", difftime(t_stop, t_start));

    free_mem(hash_dt);
    return 0;
}


static unsigned int
hash_func(char * chrom_s)
{
    if (!strcmp(chrom_s, "X"))
        return 99;
    else if (!strcmp(chrom_s, "Y"))
        return 98;
    else if (!strcmp(chrom_s, "XY"))
        return 97;
    else if (isdigit(chrom_s[0]))
        return atoi(chrom_s);
    else{
        fprintf(stderr, "error, chrom %s not recognized.\n", chrom_s);
        exit(EXIT_FAILURE);
    }
}


static const char *
unhash_func(unsigned int chrom_hash)
{
    char * chrom_s = (char *)malloc(10);

    if (chrom_hash == 99)
        chrom_s = "X";
    else if (chrom_hash == 98)
        chrom_s = "Y";
    else if (chrom_hash == 97)
        chrom_s = "XY";
    else
        sprintf(chrom_s, "%d", chrom_hash);
    return chrom_s;
    
}


static void
structure_gene_data(struct HASH_NODE * hash_dt, const char * gene_list_name, bool head, unsigned int extend_len)
{
    FILE * f_in = fopen(gene_list_name, "r");
    char line_buffer[LINE_BUFFER_LEN];
    char field_buffer[FIELD_BUFFER_LEN];
    unsigned int i = 0, j = 0, k = 0;
    char  gene_name[20];
    unsigned char chrom_hash = 0;
    unsigned long gene_start = 0;
    unsigned long gene_end = 0;
    struct GENE_NODE * new;
    struct GENE_NODE * this;
    unsigned long this_start = 0;
    unsigned long * gene_list;
    unsigned long most_right = 0;

    if (f_in == NULL){
        fprintf(stderr, "error, open gene list file failed.\n");
        exit(EXIT_FAILURE);
    }
    if (head){
        fgets(line_buffer, LINE_BUFFER_LEN, f_in);
        if (line_buffer[strlen(line_buffer) - 1] != '\n'){
            fprintf(stderr, "error, line_buffer overflow.\n");
            exit(EXIT_FAILURE);
        }
    }

    while (fgets(line_buffer, LINE_BUFFER_LEN, f_in) != NULL){
        j = 0;
        k = 0;
        for (i = 0; line_buffer[i] != '\0'; ++i){
            if (line_buffer[i] == ' ' || line_buffer[i] == '\t' || line_buffer[i] == '\n'){
                field_buffer[j] = '\0';
                j = 0;
                switch (k){
                case 0:
                    chrom_hash = hash_func(field_buffer);
                    break;
                case 1:
                    gene_start = atol(field_buffer);
                    break;
                case 2:
                    gene_end = atol(field_buffer);
                    break;
                case 3:
                    strcpy(gene_name, field_buffer);
                    break;
                default:
                    break;
                }
                ++k;
            } else{
                if (j > FIELD_BUFFER_LEN - 2){
                    fprintf(stderr, "error, field buffer overflow.\n");
                    exit(EXIT_FAILURE);
                }
                field_buffer[j] = line_buffer[i];
                ++j;
            }
        }
        if (line_buffer[i - 1] != '\n'){
            fprintf(stderr, "error2, line buffer overflow.\n");
            exit(EXIT_FAILURE);
        }

        new = (struct GENE_NODE *)malloc(sizeof(struct GENE_NODE));
        (new -> gene_pos_s_e)[0] = gene_start;
        (new -> gene_pos_s_e)[1] = gene_end;
        strcpy(new -> gene_name, gene_name);
        new -> snp_pointer = NULL;
        new -> snp_pointer_tail = NULL;
        new -> prev = NULL;
        new -> next = NULL;
        hash_dt[chrom_hash].gene_num += 1;

        if (hash_dt[chrom_hash].gene_dt == NULL){
            hash_dt[chrom_hash].gene_dt = new;
            hash_dt[chrom_hash].gene_dt_tail = new;
        }else{
            this = hash_dt[chrom_hash].gene_dt;
            while (1){
                if (this == NULL)
                    break;
                this_start = (this -> gene_pos_s_e)[0];
                if (this_start > gene_start){
                    if (this -> prev == NULL){
                        hash_dt[chrom_hash].gene_dt = new;
                        new -> next = this;
                        this -> prev = new;
                        break;
                    } else{
                        (this -> prev) -> next = new;
                        new -> prev = this -> prev;
                        new -> next = this;
                        this -> prev = new;
                        break;
                    }
                } else{
                    if (this -> next == NULL){
                        this -> next = new;
                        new -> prev = this;
                        hash_dt[chrom_hash].gene_dt_tail = new;
                        break;
                    } else{
                        this = this -> next;
                    }
                }
            }
        }
    }
    fclose(f_in);

    for (i = 0; i < HASH_LEN; ++i){
        if (hash_dt[i].gene_dt){
            gene_list = (unsigned long *)malloc(sizeof(unsigned long) * hash_dt[i].gene_num * 3);
            this = hash_dt[i].gene_dt;
            j = 0;
            most_right = 0;
            while (this){
                gene_start = ((long)((this -> gene_pos_s_e)[0] - extend_len) > 0)? (this -> gene_pos_s_e)[0] - extend_len: 1;
                gene_end = (this -> gene_pos_s_e)[1] + extend_len;
                if (gene_end > most_right)
                    most_right = gene_end;
                gene_list[j] = gene_start;
                ++j;
                gene_list[j] = gene_end;
                ++j;
                gene_list[j] = most_right;
                ++j;
                this = this -> next;
            }
            hash_dt[i].gene_order_list = gene_list;
        }
    }
    return;
}


static void
structure_snp_data(struct HASH_NODE * hash_dt, const char * snp_list_name, bool head, \
    double max_af, unsigned int min_sample)
{
    FILE * f_in = fopen(snp_list_name, "r");
    char line_buffer[LINE_BUFFER_LEN];
    char field_buffer[FIELD_BUFFER_LEN];
    unsigned int i = 0, j = 0, k =0;
    unsigned long snp_pos = 0;
    unsigned char chrom_hash = 0;
    unsigned long sample_num = 0;
    double af = 0.0;
    double p_value = 0.0;
    struct SNP_NODE * new;

    if (f_in == NULL){
        fprintf(stderr, "error, open snp list file failed.\n");
        exit(EXIT_FAILURE);
    }
    if (head){
        fgets(line_buffer, LINE_BUFFER_LEN, f_in);
        if (line_buffer[strlen(line_buffer) - 1] != '\n'){
            fprintf(stderr, "error, line_buffer overflow.\n");
            exit(EXIT_FAILURE);
        }
    }

    while (fgets(line_buffer, LINE_BUFFER_LEN, f_in) != NULL){
        j = 0;
        k = 0;
        for (i = 0; line_buffer[i] != '\0'; ++i){
            if (line_buffer[i] == ' ' || line_buffer[i] == '\t' || line_buffer[i] == '\n'){
                field_buffer[j] = '\0';
                j = 0;
                switch (k){
                case 0:
                    chrom_hash = hash_func(field_buffer);
                    break;
                case 2:
                    snp_pos = atol(field_buffer);
                    break;
                case 5:
                    sample_num = atol(field_buffer);
                    break;
                case 6:
                    af = atof(field_buffer);
                    break;
                case 12:
                    p_value = atof(field_buffer);
                    break;
                default:
                    break;
                }
                ++k;
            } else {
                if (j > FIELD_BUFFER_LEN - 2){
                    fprintf(stderr, "error, field buffer overflow.\n");
                    exit(EXIT_FAILURE);
                }
                field_buffer[j] = line_buffer[i];
                ++j;
            }
        }
        if (line_buffer[ i - 1] != '\n'){
            fprintf(stderr, "error, line buffer overflow.\n");
            exit(EXIT_FAILURE);
        }

        af = (af < 0.5)? af: 1.0 -af;
        sample_num = af * sample_num * 2;

        if ((af - max_af) <= 1e-15 && sample_num >= (unsigned long)min_sample && \
            hash_dt[chrom_hash].gene_dt != NULL && \
            judge_snp_in_gene(snp_pos, hash_dt[chrom_hash].gene_order_list, hash_dt[chrom_hash].gene_num)){
            new = (struct SNP_NODE *)malloc(sizeof(struct SNP_NODE));
            new -> next = NULL;
            new -> snp_pos = snp_pos;
            new -> weight = pow(BETA_DENSITY(af, 1, 25) / BETA_DENSITY(af, 0.5, 0.5), 2);
            new -> p_value = p_value;

            if (hash_dt[chrom_hash].snp_dt == NULL){
                hash_dt[chrom_hash].snp_dt = new;
                hash_dt[chrom_hash].snp_dt_tail = new;
            } else {
                hash_dt[chrom_hash].snp_dt_tail -> next = new;
                hash_dt[chrom_hash].snp_dt_tail = new;
            }
        }
    }
    fclose(f_in);
    return;
}


static bool
judge_snp_in_gene(unsigned long snp_pos, unsigned long * gene_order_list, unsigned long list_len)
{
    list_len = list_len * 3;
    unsigned long most_right = 0;
    unsigned long i = 0;

    for (i = 0; i < list_len; i += 3){
        if (snp_pos < gene_order_list[i]){
            if (snp_pos > most_right){
                return false;
            }else{
                return true;
            }
        }
        most_right = gene_order_list[i + 2];
    }

    if (snp_pos > most_right)
        return false;
    else
        return true;
}


static void *
thread_worker_1(void * hash_dt)
{
    struct HASH_NODE          hash_node;
    struct GENE_NODE        * gene_dt_head;
    struct GENE_NODE        * gene_dt_pt;
    struct GENE_NODE        * gene_dt_tail;
    struct GENE_NODE        * gene_dt_backtrack;
    struct SNP_POINTER_NODE * snp_pt_node;
    struct SNP_NODE         * snp_dt_head;
    unsigned long           * gene_order_list;
    unsigned long             list_len;
    unsigned long             i = 0;
    unsigned long             gene_index = 0;
    unsigned long             snp_pos = 0;
    int                       out_mark  = 0;

    hash_node = *((struct HASH_NODE *)hash_dt);
    gene_dt_head = hash_node.gene_dt;
    gene_order_list = hash_node.gene_order_list;
    list_len = hash_node.gene_num * 3;
    snp_dt_head = hash_node.snp_dt;
    gene_dt_tail = hash_node.gene_dt_tail;

    while(snp_dt_head){
        snp_pos = snp_dt_head -> snp_pos;
        out_mark = 1;
        gene_dt_pt = gene_dt_head;
        gene_index = 0;

        for (i = 0; i < list_len; i += 3){
            if (snp_pos < gene_order_list[i]){
                out_mark = 0;
                gene_dt_backtrack = gene_dt_pt -> prev;
                for (; gene_index > 0; gene_index -= 3){
                    if (snp_pos <= gene_order_list[gene_index - 1]){
                        if(snp_pos <= gene_order_list[gene_index -2]){
                            snp_pt_node = (struct SNP_POINTER_NODE *)malloc(sizeof(struct SNP_POINTER_NODE));
                            snp_pt_node -> snp_pt = snp_dt_head;
                            snp_pt_node -> next = NULL;
                            if (gene_dt_backtrack -> snp_pointer == NULL){
                                gene_dt_backtrack -> snp_pointer = snp_pt_node;
                                gene_dt_backtrack -> snp_pointer_tail = snp_pt_node;
                            } else{
                                (gene_dt_backtrack -> snp_pointer_tail) -> next = snp_pt_node;
                                gene_dt_backtrack -> snp_pointer_tail = snp_pt_node;
                            }
                        }
                        gene_dt_backtrack = gene_dt_backtrack -> prev;
                    } else{
                        break;
                    }
                }
                break;
            }
            gene_dt_pt = gene_dt_pt -> next;
            gene_index += 3;
        }

        if (out_mark){
            gene_dt_backtrack = gene_dt_tail;
            for (gene_index = list_len; gene_index > 0; gene_index -= 3){
                if (snp_pos <= gene_order_list[gene_index - 1]){
                    if (snp_pos <= gene_order_list[gene_index - 2]){
                        snp_pt_node = (struct SNP_POINTER_NODE *)malloc(sizeof(struct SNP_POINTER_NODE));
                        snp_pt_node -> snp_pt = snp_dt_head;
                        snp_pt_node -> next = NULL;
                        if (gene_dt_backtrack -> snp_pointer == NULL){
                            gene_dt_backtrack -> snp_pointer = snp_pt_node;
                            gene_dt_backtrack -> snp_pointer_tail = snp_pt_node;
                        } else{
                            (gene_dt_backtrack -> snp_pointer_tail) -> next = snp_pt_node;
                            gene_dt_backtrack -> snp_pointer_tail = snp_pt_node;
                        }
                    }
                }else{
                    break;
                }
                gene_dt_backtrack = gene_dt_backtrack -> prev;
            }
        }
        snp_dt_head = snp_dt_head -> next;
    }
    return NULL;
}


static void
print_res_file(struct HASH_NODE * hash_dt, const char * f_out_name)
{
    FILE * f_out = fopen(f_out_name, "w");
    unsigned char i = 0;
    struct GENE_NODE * gene_dt;
    struct SNP_POINTER_NODE * snp_pt_node;
    double weight_sum = 0.0;
    double weight_ave = 0.0;
    double c_value_sum = 0.0;
    double p_value;
    double tmp;
    double cauchy = 0.0;
    unsigned long snp_amount;

    fprintf(f_out, "CHR\tGENE\tSTART\tEND\tSNP_NUM\tP_ACAT\n");

    for (i = 0; i < HASH_LEN; i++){
        if(hash_dt[i].gene_dt != NULL && hash_dt[i].snp_dt != NULL){
            gene_dt = hash_dt[i].gene_dt;
            while(1){
                if (gene_dt == NULL)
                    break;
                if (gene_dt -> snp_pointer != NULL){
                    weight_sum = 0.0;
                    c_value_sum = 0.0;
                    snp_amount = 0;

                    snp_pt_node = gene_dt -> snp_pointer;
                    while(1){
                        if (snp_pt_node == NULL)
                            break;
                        weight_sum += (snp_pt_node -> snp_pt) -> weight;
                        snp_amount++;
                        snp_pt_node = snp_pt_node -> next;
                    }

                    snp_pt_node = gene_dt -> snp_pointer;
                    while(1){
                        if (snp_pt_node == NULL)
                            break;
                        p_value = (snp_pt_node -> snp_pt) -> p_value;
                        weight_ave = ((snp_pt_node -> snp_pt) -> weight) / weight_sum;
                        if (p_value >= 1e-15){
                            tmp = (0.5 - p_value) * M_PI;
                            c_value_sum += weight_ave * (sin(tmp) / cos(tmp));
                        }else{
                            c_value_sum += (weight_ave / p_value) * M_PI;
                        }
                        snp_pt_node = snp_pt_node -> next;
                    }
                    if (fabs(c_value_sum) > 1){
                        cauchy = PCAUCHY2(c_value_sum);
                        cauchy = (c_value_sum > 0)? cauchy: 1 + cauchy;

                    } else{
                        cauchy  = 1.0L - PCAUCHY(c_value_sum);
                    }

                    fprintf(f_out, "%s\t%s\t%ld\t%ld\t%ld\t%le\n", unhash_func(i), gene_dt -> gene_name, \
                        (gene_dt -> gene_pos_s_e)[0], (gene_dt -> gene_pos_s_e)[1],  snp_amount, cauchy);
                }
                gene_dt = gene_dt -> next;
            }
        }
    }
    fclose(f_out);
    return;
}


static void
free_mem(struct HASH_NODE * hash_dt)
{
    unsigned char i;
    struct GENE_NODE * gene_dt;
    struct GENE_NODE * gene_dt_tmp;
    struct SNP_POINTER_NODE * snp_pt_node;
    struct SNP_POINTER_NODE * snp_pt_tmp;
    struct SNP_NODE * snp_dt;
    struct SNP_NODE * snp_dt_tmp;

    for (i = 0; i < HASH_LEN; i++){
        if (hash_dt[i].gene_dt){
            hash_dt[i].gene_dt_tail = NULL;
            free(hash_dt[i].gene_order_list);
            gene_dt = hash_dt[i].gene_dt;
            while(gene_dt){
                if (gene_dt -> snp_pointer){
                    gene_dt -> snp_pointer_tail = NULL; 
                    snp_pt_node = gene_dt -> snp_pointer;
                    while(snp_pt_node){
                        snp_pt_tmp = snp_pt_node -> next;
                        free(snp_pt_node);
                        snp_pt_node = snp_pt_tmp;
                    }
                }
                gene_dt_tmp = gene_dt -> next;
                if (gene_dt_tmp != NULL)
                    gene_dt_tmp -> prev = NULL;
                free(gene_dt);
                gene_dt = gene_dt_tmp;
            }
        }
        if (hash_dt[i].snp_dt){
            hash_dt[i].snp_dt_tail = NULL;
            snp_dt = hash_dt[i].snp_dt;
            while(snp_dt){
                snp_dt_tmp = snp_dt -> next;
                free(snp_dt);
                snp_dt = snp_dt_tmp;
            }
        }
    }
    return;
}


