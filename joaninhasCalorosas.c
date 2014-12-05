/*
 ====================
 JoaninhasCalorosas.c
 ====================

 ==============================================================================
 MAC0431 -- 28/11/2014 -- IME/USP, -- Prof. Marco Dimas Gubitoso
 Author      : Marcello Souza de Oliveira, 6432692
			 : Francisco Alexandre da Silva Gomes dos Santos, 794650
			 : Mauricio Santana, 7991170

 Course      : Bacharelado em Ciencias da Computacao

 Name        : JoaninhasCalorosas.c
 Version     :
 Copyright   :
 Description : Joaninhas Calorosas: simulacao...
   ...

 Compiler         : gcc linux 4.6.3
 Compiler Options : -Wall -ansi -pedantic -O2 -U_FORTIFY_SOURCE
 Link Options     : -fopenmp -lm
 Editor           : Sublime Text 2, Geany 0.21, Code-Blocks 13.12;
 S.O.             : Linux
 ==============================================================================

 ==============================================================================
 References:
 [1] http://www.redblobgames.com/grids/hexagons/
     com “even-r” horizontal layout.
 ==============================================================================

 ==============================================================================
 How to compile:
 $ gcc -o joaninhasCalorosas joaninhasCalorosas.c -fopenmp -lm
 $ ./joaninhasCalorosas <params>

 $ rm joaninhasCalorosas ; make -k joaninhasCalorosas && /usr/bin/time ./joaninhasCalorosas 100 100 30 0 22.765 18 27 0.1 30 0.1 40 10000 2

 ==============================================================================
*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h> /* 'linking' com '-lm' flag: gcc fontes.c bibs.h -o exec -lm */
#include <time.h>
/*#include <sys/time.h>*/
#include <unistd.h>

/*#include <mpi.h>*/
#include <omp.h>

/********
 * MACROS
 *******/

#define TOL 1.00

#ifndef M_PI
	#define M_PI 3.14159265358979323846  /* pi in math.h */
#endif

#define min( X, Y) (((X) < (Y)) ? (X) : (Y))

#define max( a, b) \
	({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; })

/**********
 * DEBUGING
 *********/

/*
 * Habilita o modo de depuração.
 */
//#define DEBUG

/*
 * Habilita a impressao de resultados no arquivo jc-out.txt.
 */
#define PRINTBUGS

/*
 * Habilita a impressao do tempo decorrido entre o inicio e fim do trecho de
 * simulacao (sem escrita no arq. saida).
 */
#define PRINTTIME

#define NOCALC -1

#ifdef DEBUG
	#define MSG_DBG( fmt, ...) fprintf( stderr, "\t[MSG_DBG:] %s (%s:%d): " \
										fmt"\n", __FUNCTION__, __FILE__, \
										__LINE__, ##__VA_ARGS__)
#else
   #define MSG_DBG( fmt, ...) {;}
#endif

/*
 * Ref.:
 * http://cboard.cprogramming.com/linux-programming/113208-portable-method-
 * accessing-cpu-time-given-process.html
 * TELLTIME;
 */
#ifdef _POSIX_CPUTIME
	#define STARTTIME() diff = clock()
	#define TELLTIME() diff = clock() -diff; printf( "\n\tTime Elapsed = %f\n\n", (double) ((double) diff/CLOCKS_PER_SEC))
	#define TELLTIME2() fprintf( stdout, "Time Elapsed: %G seconds.\t(%s: %d)\n", \
							  (double) ((double) clock( )/CLOCKS_PER_SEC), \
							  __FUNCTION__, __LINE__);
#endif

/*********
 * STRUCTS
 ********/

typedef enum { VACANT = 0, LADYBUG, HOT, COLD } occupation;
typedef enum { D_HOT = 0, D_COLD = 1 } disType;
typedef enum { EVEN, ODD } parity;

/*
 * index: possibilita testar ao final dos calculos se duas joaninhas desejam
 * se mover para o mesmo hex (devera prevalecer a com maior diferenca de
 * energia para com a energia do hex destino, conforme enunciado), sem a
 * necessidade de percorrer todo o vetor de joaninhas (mas pela indexação).
 * index: indice para o vetor de joaninhas
 * sem: semente atual
 */
typedef struct hex {
	int index;
	occupation elem;
	int sem;
} hex;
typedef hex **hexgrid;

typedef struct position {
	int row;
	int col;
} position;

typedef struct dposition {
	position pos;
	double energy;
	unsigned int bornscycle;
} dposition;

typedef struct eposition {
	position pos;
	double energy;
} eposition;

/*
 * amount: 0 - 6 vizinhos possiveis (considera-se apenas hex's VAZIA's).
 */
typedef struct bug {
	eposition orig;
	eposition dest;
	eposition vacneighbor[6];		/* vacant neighbors array */
	int amount;
	int movesflag;
} bug;

typedef struct entData {
	unsigned int A, L, j, s, nc, nf, T, P;
	float pc, pf,					/* probs em decimal */
		  th_min, th_max;
	double C;
} entData;

struct disturbs {
	dposition *pos;
	int begin;
	int end;
	int size;
	int ncycle;
	float prob;
};
typedef struct disturbs disturbs;

/******************
 * GLOBAL VARIABLES
 *****************/

//static unsigned int cycle;

clock_t diff;

/************
 * PROTOTYPES
 ***********/

void *mallocX( unsigned int nbytes);

hex **createHexGrid( const entData *data);

void destroyHexGrid( hex *H[], const entData *data);

bug *createBugs( hex *H[], entData *data);

void showBugs( bug ladybug[], const entData *data);

void killBugs( bug *ladybug);

void QUEUEinit( disturbs dis[], const entData *data);

void QUEUEput( disturbs dis[], disType elem, int r_row, int q_col, unsigned int bornscycle);

void QUEUEupdate( hex *H[], disturbs dis[], unsigned int cycle);

void QUEUEfree( disturbs dis[]);

void generateDisturbs( hex *H[], disturbs dis[], entData *data, unsigned int cycle);

void neighbornsBugs( hex *H[], int L, int A, bug ladybug[], int nbugs);

void bugsMove( hex *H[], const entData *data, bug ladybug[]);

double partialEnergy( position p, position q, occupation q_elem, double C);

void disturbsEnergy( eposition *orig, disturbs dis[], double C);

void calcEnergy( hex *H[], disturbs dis[], bug ladybug[], const entData *data);

void nextCycle( hex *H[], disturbs dis[], entData *data, unsigned int cycle);

int startSimulation( hex *H[], bug ladybug[], entData *data);

/*
 ==============================================================================
*/

/***************
 * MAIN FUNCTION
 **************/

int main( int argc, char *argv[]) {
	entData data;
	bug *ladybug;			/* Bug's positions vector */
	hex **H;				/* Hex Matrix (Hex Grid) */

	if (argc == 1) {
	/*
		data.L = 500;
		data.A = 500;
		data.j = 10;
		data.s = 0;
		data.C = 1;
		data.th_min = 1.25;
		data.th_max = 2.0;
		data.pc = 0.4;
		data.nc = 2;
		data.pf = 0.3;
		data.nf = 1;
		data.T = 20;
		data.P = 1;
	*/
		data.L = 100;
		data.A = 100;
		data.j = 100;
		data.s = 0;
		data.C = 22.765;
		data.th_min = 18;
		data.th_max = 27;
		data.pc = 0.1;
		data.nc = 30;
		data.pf = 0.1;
		data.nf = 40;
		data.T = 10000;
		data.P = 2;
	}
	else if (argc < 14) {
		fprintf( stderr,
				 "%d - Modo de uso:\n\t$ .%s L A j s C th_min th_max pc nc pf nf T P\n\n",
				 -1, argv[0]);
		fprintf( stderr, "onde:\n"
				 " L : largura do hex grid;\n"
				 " A : altura do hex grid;\n"
				 " j : num. de joaninhas;\n"
				 " s : semente para gerador aleatorio;\n"
				 " C : constante de emissao de calor de uma joaninha;\n");
		fprintf( stderr,
				 " th_min : menor temperatura que a joaninha considera confortavel;\n"
				 " th_max : maior temperatura que a joaninha considera confortavel;\n"
				 " pc : prob./hex de aparecer uma fonte de calor;\n"
				 " nc : duracao em ciclos da fonte de calor;\n"
				 " pf : prob./hex de aparecer uma fonte de frio;\n"
				 " nf : duracao em ciclos da fonte de calor;\n"
				 " T : num. de ciclos da simulacao;\n"
				 " P : num. de threads na simulacao.\n\n");
		exit( EXIT_FAILURE);
	}
	else {
/*
		scanf( " %u %u %u %u %lf %lf %lf %f %u %f %u %u %u",
			   &data.L, &data.A, &data.j, &data.s, &data.C, &data.th_min,
			   &data.th_max, &data.pc, &data.nc, &data.pf, &data.nf,
			   &data.T, &data.P);
*/
		data.L = atoi( argv[1]);
		data.A = atoi( argv[2]);
		data.j = atoi( argv[3]);
		data.s = atoi( argv[4]);
		data.C = atof( argv[5]);
		data.th_min = atof( argv[6]);
		data.th_max = atof( argv[7]);
		data.pc = atof( argv[8]);
		data.nc = atoi( argv[9]);
		data.pf = atof( argv[10]);
		data.nf = atoi( argv[11]);
		data.T = atoi( argv[12]);
		data.P = atoi( argv[13]);
	}

	if (data.th_min > data.th_max)
		data.th_max = data.th_min;

	omp_set_num_threads( data.P);

	H = createHexGrid( &data);
	ladybug = createBugs( H, &data);

	#ifdef PRINTTIME
		STARTTIME();
	#endif
	startSimulation( H, ladybug, &data);
	#ifdef PRINTTIME
		TELLTIME();
	#endif

	#ifdef PRINTBUGS
		showBugs( ladybug, &data);
	#endif

	killBugs( ladybug);
	destroyHexGrid( H, &data);

	return EXIT_SUCCESS;
}

/***********
 * FUNCTIONS
 **********/

/*
 * Versão da func. malloc com tratamento de erro.
 */
void *mallocX( unsigned int nbytes) {
	void *ptr;
	ptr = malloc( nbytes);
	if (ptr == NULL) {
		fprintf( stderr, "%d -%d- A memoria estourou! Socorro! "
				 "malloc devolveu NULL! (%d bytes)\n", -1, __LINE__, nbytes);
		exit( EXIT_FAILURE);
	}
	return ptr;
}

/*
 * Matriz Hexagonal. A distância ("euclidiana") entre os centros de 2
 * hexágonos vizinhos é de uma unidade.
 * “even-r” horizontal layout.
 */
hex **createHexGrid( const entData *data) {
	int a, l;
	hex **H = (hex **) malloc( data->A * sizeof (struct hex *));
	if (H == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
				 "malloc devolveu NULL! (%d bytes)\n",
				 -1, (int) (data->A * sizeof (struct hex *)));
		exit( EXIT_FAILURE);
	}
	#pragma omp parallel
	{
		#pragma omp for private (a, l)
		for (a = 0; a < data->A; ++a) {
			H[a] = (hex *) malloc( data->L * sizeof (struct hex));
			if (H[a] == NULL) {
				fprintf( stderr, "%d -%d- A memoria estourou! Socorro! "
						 "malloc devolveu NULL! (%d bytes)\n",
						 -1, __LINE__, (int) (data->L * sizeof (struct hex)));
				exit( EXIT_FAILURE);
			}
			for (l = 0; l < data->L; ++l) {
				H[a][l].elem = VACANT;
				H[a][l].index = data->j;
				H[a][l].sem = ((a + 1) * data->s + l) % RAND_MAX;
			}
		}
	}
	return H;
}

void destroyHexGrid( hex *H[], const entData *data) {
	int a;
	if (H != NULL) {
		for (a = 0; a < data->A; ++a) {
			free (H[a]);
		}
		free (H);
	}
	return;
}

bug *createBugs( hex *H[], entData *data) {
	int i, k, r, q;
	bug *p = (bug *) malloc( data->j * sizeof (struct bug));
	if (p == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
				 "malloc devolveu NULL! (%d bytes)\n",
				 -1, (int) (data->j * sizeof (struct bug)));
		exit( EXIT_FAILURE);
	}
	for (i = 0, k = 0; i < data->j;) {
		r = (int) (((double) rand_r( &data->s) / (RAND_MAX + 1.0)) * data->A);
		q = (int) (((double) rand_r( &data->s) / (RAND_MAX + 1.0)) * data->L);
		if (H[r][q].elem == VACANT) {
			H[r][q].elem = LADYBUG;
			H[r][q].index = k;
			p[k].dest.pos.row = p[k].orig.pos.row = r;
			p[k].dest.pos.col = p[k].orig.pos.col = q;
			p[k].orig.energy = 0.0;
			p[k].amount = NOCALC;
			p[k].movesflag = 0;
			i++; k++;
		}
	}
	return p;
}

void showBugs( bug ladybug[], const entData *data) {
	int b;
	FILE *out = fopen( "jc-out.txt", "w");
	/*FILE *out = stdout;*/
	if (out != NULL) {
		/*fprintf( out, "Joaninha J( row, col) E energy\n\n");*/
		for (b = 0; b < data->j; ++b) {
			fprintf( out, "  %5d    %5d \t %.10f\n", ladybug[b].orig.pos.row, ladybug[b].orig.pos.col, ladybug[b].orig.energy);
		}
		fprintf( out, "\n");
		fclose( out);
	}
}

void killBugs( bug *ladybug) {
	if (ladybug != NULL)
		free (ladybug);
	return;
}

/*
 ==============================================================================
*/

void QUEUEinit( disturbs dis[], const entData *data) {
	int I;
	if (dis == NULL) return;
	dis[D_HOT].prob = data->pc;
	dis[D_HOT].ncycle = data->nc;
	dis[D_COLD].prob = data->pf;
	dis[D_COLD].ncycle = data->nf;
//	#pragma omp parallel
//	{
//		#pragma omp for private (I)
		for (I = D_HOT; I <= D_COLD; ++I) {
			dis[I].size = ((data->L * data->A) - data->j +1);
			dis[I].begin = 0;
			dis[I].end = 0;
			dis[I].pos = (dposition *) malloc( dis[I].size * sizeof (struct dposition));
			if (dis[I].pos == NULL) {
				fprintf( stderr, "%d -%d- A memoria estourou! Socorro! "
						 "malloc devolveu NULL! (%d bytes)\n",
						 -1, __LINE__, (int) (dis[I].size * sizeof (struct dposition)));
				exit( EXIT_FAILURE);
			}
		}
//	}
}

/*
 *
 */
void QUEUEput( disturbs dis[], disType elem, int r_row, int q_col, unsigned int bornscycle) {
	int k;
	#pragma omp critical
	k = dis[elem].end;
	dis[elem].end = (dis[elem].end + 1) % dis[elem].size;
	#pragma omp parallel
	/* verifica se dis[elem] esta cheio e se estiver, aloca novas posicoes */
	if (dis[elem].end == dis[elem].begin) {
		printf( "\n%d - Socorro! Fila vai transbordar!\n", __LINE__);
		exit( EXIT_FAILURE);
	}
	dis[elem].pos[k].pos.row = r_row;
	dis[elem].pos[k].pos.col = q_col;
	dis[elem].pos[k].energy = 0.0;
	dis[elem].pos[k].bornscycle = bornscycle;
	/*(dis[elem].end + 1 == dis[elem].size) ? dis[elem].end = 0 : dis[elem].end++;*/
}

void QUEUEupdate( hex *H[], disturbs dis[], unsigned int cycle) {
	int I;
	for (I = D_HOT; I < D_COLD; ++I) {
		if (dis[I].begin != dis[I].end)
			while (dis[I].pos[dis[I].begin].bornscycle + dis[I].ncycle <= cycle && dis[I].begin != dis[I].end) {
				H[dis[I].pos[dis[I].begin].pos.row][dis[I].pos[dis[I].begin].pos.col].elem = VACANT;
				(dis[I].begin + 1 == dis[I].size) ? dis[I].begin = 0 : dis[I].begin++;
			}
	}
}

void QUEUEfree( disturbs dis[]) {
	int I;
	for (I = D_HOT; I < D_COLD; ++I){
		if (dis[I].pos != NULL)
			free( dis[I].pos);
	}
}

/*
 ==============================================================================
*/

void generateDisturbs( hex *H[], disturbs dis[], entData *data, unsigned int cycle) {
	int a, l, k;
	#pragma omp parallel shared (H, dis, k)
	{
		#pragma omp for private (a, l, k)
		for (a = 0; a < data->A; ++a) {
			for (l = 0; l < data->L; ++l) {
				/*H[a][l].sem = ((a + 1) * H[a][l].sem + l) % RAND_MAX;*/
				if (H[a][l].elem == VACANT) {
					if (rand_r( &H[a][l].sem) / (RAND_MAX + 1.0) <= data->pc) {
						H[a][l].elem = HOT;
						#pragma omp critical
						k = dis[D_HOT].end;
						dis[D_HOT].end = (dis[D_HOT].end + 1) % dis[D_HOT].size;
						#pragma omp parallel
						/* verifica se dis[D_HOT] esta cheio e se estiver, aloca novas posicoes */
						if (dis[D_HOT].end == dis[D_HOT].begin) {
							printf( "\n%d - Socorro! Fila vai transbordar!\n", __LINE__);
							exit( EXIT_FAILURE);
						}
						dis[D_HOT].pos[k].pos.row = a;
						dis[D_HOT].pos[k].pos.col = l;
						dis[D_HOT].pos[k].energy = 0.0;
						dis[D_HOT].pos[k].bornscycle = cycle;
					}
					else if (rand_r( &H[a][l].sem) / (RAND_MAX + 1.0) <= data->pf) {
						H[a][l].elem = COLD;
						#pragma omp critical
						k = dis[D_COLD].end;
						dis[D_COLD].end = (dis[D_COLD].end + 1) % dis[D_COLD].size;
						#pragma omp parallel
						/* verifica se dis[D_COLD] esta cheio e se estiver, aloca novas posicoes */
						if (dis[D_COLD].end == dis[D_COLD].begin) {
							printf( "\n%d - Socorro! Fila vai transbordar!\n", __LINE__);
							exit( EXIT_FAILURE);
						}
						dis[D_COLD].pos[k].pos.row = a;
						dis[D_COLD].pos[k].pos.col = l;
						dis[D_COLD].pos[k].energy = 0.0;
						dis[D_COLD].pos[k].bornscycle = cycle;
					}
				}
			}
		}
	}
	return;
}

/*
 * Calcula vizinhos vazios de cada joaninha.
 *
 * Even-r layout:
 * paritying=0 => even-r (i=0,2,4,...):
 * 		E=r,q+1;NE=r-1,q+1;NW=r-1,q;W=r,q-1;SW=r+1,q;SE=r+1,q+1
 * paritying=1 => odd-r  (i=1,3,5,...):
 * 		E=r,q+1;NE=r-1,q;NW=r-1,q-1;W=r,q-1;SW=r+1,q-1;SE=r+1,q
 */
void neighbornsBugs( hex *H[], int L, int A, bug ladybug[], int nbugs) {
	int i, j, k,
		paritying,
		row, col;
	int directions[2][6][2] = {
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } },
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } }
	}; /* even-r layout - E,NE,NW,W,SW,SE */
	for (j = 0; j < nbugs; ++j) {
		paritying = ladybug[j].orig.pos.row & 1;
		for (i = 0, k = 0; i < 6; ++i) {
			row = ladybug[j].orig.pos.row + directions[paritying][i][0];
			col = ladybug[j].orig.pos.col + directions[paritying][i][1];
			if (row >= 0 && row < A && col >= 0 && col < L) {
				if (H[row][col].elem == VACANT) {
					ladybug[j].vacneighbor[k].pos.row = row;
					ladybug[j].vacneighbor[k].pos.col = col;
					ladybug[j].vacneighbor[k].energy = 0.0;
					ladybug[j].amount = ++k;
				}
			}
		}
	}
}

void bugsMove( hex *H[], const entData *data, bug ladybug[]) {
	int i, b1, b2,
		paritying,
		row, col;
	int directions[2][6][2] = {
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } },
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } }
	}; /* even-r layout - E,NE,NW,W,SW,SE */
//	#pragma omp parallel shared (H, ladybug)
//	{
//		#pragma omp for private (b1, i, row, col, paritying)
		for (b1 = 0; b1 < data->j; ++b1) {
			if (ladybug[b1].movesflag == 1) {
				for (i = 0; i < 6; ++i) {
					paritying = ladybug[b1].dest.pos.row & 1;
					row = ladybug[b1].dest.pos.row + directions[paritying][i][0];
					col = ladybug[b1].dest.pos.col + directions[paritying][i][1];
					if (row >= 0 && row < data->A && col >= 0 && col < data->L &&
					  H[row][col].elem == LADYBUG && ladybug[H[row][col].index].movesflag == 1 && H[row][col].index != b1) {
						if (ladybug[H[row][col].index].dest.pos.row == ladybug[b1].dest.pos.row &&
						  ladybug[H[row][col].index].dest.pos.col == ladybug[b1].dest.pos.col) {
							if (abs( ladybug[H[row][col].index].orig.energy - ladybug[b1].dest.energy)
							  > abs( ladybug[b1].orig.energy - ladybug[b1].dest.energy)) {
								ladybug[b1].movesflag = 0;
							}
							else if (abs( ladybug[H[row][col].index].orig.energy - ladybug[b1].dest.energy)
							  < abs( ladybug[b1].orig.energy - ladybug[b1].dest.energy)) {
								ladybug[H[row][col].index].movesflag = 0;
							}
							else {
								ladybug[b1].movesflag = 0;
								ladybug[H[row][col].index].movesflag = 0;
							}
						}
					}
				}
			}
		}
//		#pragma omp barrier
//	}
	#pragma omp parallel shared (H, ladybug)
	{
		#pragma omp for private (b1, b2)
		for (b1 = 0; b1 < data->j; ++b1) {
			for (b2 = 0; b2 < ladybug[b1].amount; ++b2) {
				/* posso estar atualizando um vizinho ja atualizado (se alguma
				   joaninha ja se moveu pra ca). Nao ha problemas pois se foi
				   atualizada, foi com 0.0 mesmo */
				ladybug[b1].vacneighbor[b2].energy = 0.0;
			}
			ladybug[b1].amount = NOCALC;
			ladybug[b1].orig.energy = 0.0;
			ladybug[b1].dest.energy = 0.0;
			if (ladybug[b1].movesflag == 1) {
				H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].elem = VACANT;
				H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].index = data->j;
				ladybug[b1].orig.pos.row = ladybug[b1].dest.pos.row;
				ladybug[b1].orig.pos.col = ladybug[b1].dest.pos.col;
				H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].elem = LADYBUG;
				H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].index = b1;
				ladybug[b1].movesflag = 0;
			}
		}
	}
}

/*
 * Distância Euclidiana em uma Hex Grid Matrix, even-r layout:
 *      dy = y2 - y1;
 *      dx = x2 - x1;
 *      if ((int) dx & 1) {
 *              dy += ((x1 & 1) ? 0.5 : -0.5);
 *      }
 *      dist = sqrt( 0.75*dx*dx + dy*dy);
 */
double partialEnergy( position p, position q, occupation q_elem, double C) {
	double dx, dy;
	if (q_elem == VACANT)
		return 0.0;
	else if (p.row == q.row && p.col == q.col) {
		printf( "--%d--AQUI---- p(%d,%d) -> q(%d,%d) : q_elem=%d\n", __LINE__, p.row, p.col, q.row, q.col, q_elem);
		sleep( 0.5);
		return 0.0;
	}
	else {
		dx = (double) (q.col - p.col);
		dy = (double) (q.row - p.row);
		if ((int) dy & 1) {
			dx += (double) ((p.row & 1) ? 0.5 : -0.5);
		}
	}
	/* JOANINHA ou CALOR */
	return (q_elem != COLD) ? C/(0.75*dy*dy + dx*dx) : -C/(0.75*dy*dy + dx*dx);
}

/*
 * se energy != 0.0, suponho que seja um hex cuja energia ja foi calculada, e
 * nao recalculo; c.c., mesmo que ja tenha sido calculada
 * (fontes que se compensem), (re)calculo.
 * Importante: sempre que mover uma joaninha, deve-se zerar a energia de sua
 * posicao, de todos os seus vizinhos e ainda o amount (numero de vizinhos).
 */
void disturbsEnergy( eposition *orig, disturbs dis[], double C) {
	int c, I, i, k;
	double E;
	if (orig->energy == 0.0) {
		for (I = D_HOT; I < D_COLD; ++I) {
			#pragma omp parallel private (c, k)
			{
				k = ((dis[I].end - dis[I].begin) + dis[I].size) % dis[I].size;
				E = 0.0;
				#pragma omp for private (i) reduction (+: E)
				for (i = 0; i < k; ++i) {
					c = (dis[I].begin + i) % dis[I].size;
					E += partialEnergy( orig->pos, dis[I].pos[c].pos, HOT, C);
				}
			}
			orig->energy += E;
		}
	}
}

void calcEnergy( hex *H[], disturbs dis[], bug ladybug[], const entData *data) {
	int b1, b2, k, minIndex, maxIndex;
	double minE, maxE, E1, E2;
	neighbornsBugs( H, data->L, data->A, ladybug, data->j);
	#pragma omp parallel shared (H, ladybug, dis)
	{
		#pragma omp for private (b1, b2, k, minE, maxE, \
								 minIndex, maxIndex) reduction (+: E1, E2)
		for (b1 = 0; b1 < data->j; ++b1) {
			disturbsEnergy( &ladybug[b1].orig, dis, data->C);
			E1 = 0.0;
			for (b2 = 0; b2 < data->j; ++b2) {
				if (b2 != b1) {
if (ladybug[b1].orig.pos.row == ladybug[b2].orig.pos.row && ladybug[b1].orig.pos.col == ladybug[b2].orig.pos.col)
printf( "--%d--AQUI---- p(%d,%d) -> q(%d,%d)\n", __LINE__, ladybug[b1].orig.pos.row, ladybug[b1].orig.pos.col, ladybug[b2].orig.pos.row, ladybug[b2].orig.pos.col);

					E1 += partialEnergy( ladybug[b1].orig.pos, ladybug[b2].orig.pos, LADYBUG, data->C);
				}
			}
			ladybug[b1].orig.energy += E1;
			if (ladybug[b1].amount > 0) {
				disturbsEnergy( &ladybug[b1].vacneighbor[0], dis, data->C); /* atualizou: ladybug[b1].vacneighbor[0].energy */
				minE = maxE = ladybug[b1].vacneighbor[0].energy;
				minIndex = maxIndex = 0;
				for (k = 0; k < ladybug[b1].amount; ++k) {
					disturbsEnergy( &ladybug[b1].vacneighbor[k], dis, data->C); /* atualizou: ladybug[b1].vacneighbor[k].energy */
					E2 = 0.0;
					for (b2 = 0; b2 < data->j; ++b2) {
if (ladybug[b1].vacneighbor[k].pos.row == ladybug[b2].orig.pos.row && ladybug[b1].vacneighbor[k].pos.col == ladybug[b2].orig.pos.col) {
printf( ">>%d>>AQUI---- p(%d,%d) -> q(%d,%d)\n", __LINE__, ladybug[b1].vacneighbor[k].pos.row, ladybug[b1].vacneighbor[k].pos.col, ladybug[b2].orig.pos.row, ladybug[b2].orig.pos.col);
sleep(1);
}

						E2 += partialEnergy( ladybug[b1].vacneighbor[k].pos, ladybug[b2].orig.pos, LADYBUG, data->C);
					}
					ladybug[b1].vacneighbor[k].energy += E2;
					if (ladybug[b1].vacneighbor[k].energy < minE) {
						minE = ladybug[b1].vacneighbor[k].energy;
						minIndex = k;
					}
					else if (ladybug[b1].vacneighbor[k].energy > maxE) {
						maxE = ladybug[b1].vacneighbor[k].energy;
						maxIndex = k;
					}
				}
				for (k = 0; k < ladybug[b1].amount; ++k) {
					if (minE == ladybug[b1].vacneighbor[k].energy && minIndex != k) {
						/* Existe mais de um vizinho com temperatura minima
						 * --> empace caso deva procurar frio --> fico parado
						 * atualizo o index para um que nao posso usar.
						 */
						minIndex = ladybug[b1].amount;
					}
					if (maxE == ladybug[b1].vacneighbor[k].energy && maxIndex != k) {
						/* Existe mais de um vizinho com temperatura maxima
						 * --> empace caso deva procurar calor --> fico parado
						 * atualizo o index para um que nao posso usar.
						 */
						maxIndex = ladybug[b1].amount;
					}
				}
				if (ladybug[b1].orig.energy < data->th_min && maxIndex != ladybug[b1].amount) { /* procura fogo */
					ladybug[b1].dest.pos.row = ladybug[b1].vacneighbor[maxIndex].pos.row;
					ladybug[b1].dest.pos.col = ladybug[b1].vacneighbor[maxIndex].pos.col;
					ladybug[b1].dest.energy = maxE;
					ladybug[b1].movesflag = 1;
				}
				else if (ladybug[b1].orig.energy > data->th_max && minIndex != ladybug[b1].amount) { /* procura gelo */
					ladybug[b1].dest.pos.row = ladybug[b1].vacneighbor[minIndex].pos.row;
					ladybug[b1].dest.pos.col = ladybug[b1].vacneighbor[minIndex].pos.col;
					ladybug[b1].dest.energy = minE;
					ladybug[b1].movesflag = 1;
				}
			}
		}
	}
}

/*
 ==============================================================================
*/

int startSimulation( hex *H[], bug ladybug[], entData *data) {
	int cycle = 0;
	disturbs dis[2];

	QUEUEinit( dis, data);

	for (cycle = 1; cycle < data->T; ++cycle) {
//printf( "\n@--- cycle = %u ---@\n", cycle);

		generateDisturbs( H, dis, data, cycle);
		calcEnergy( H, dis, ladybug, data);
		bugsMove( H, data, ladybug);
		QUEUEupdate( H, dis, cycle);
	}
//printf( "\n@--- cycle = %u ---@\n", cycle);
	generateDisturbs( H, dis, data, cycle);
	calcEnergy( H, dis, ladybug, data);

	return EXIT_SUCCESS;
}

/*
 ==============================================================================
*/

//	#ifdef PRINT
//		printf( " %19.16G  %19.16G  %19.16G\n",
//				t_kless[0], ab->pfunc_y( t_kless[0]), y_kless[0]);*/
//	#endif

//MSG_DBG( "t_0=%f\ty_0=%f\th=%f\tt_n=%f\n", ab->t_0, ab->y_0, h, ab->t_n);

/*
int main( int argc, char *argv[]) {
	//puts( "!!!Hello World!!!"); *** prints !!!Hello World!!!

	int rank, size;
	MPI_Init( &argc, &argv);
	MPI_Comm_size( MPI_COMM_WORLD, &size);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	if (rank == 0)
		printf( "Ola, mundo. Eu sou o processo no. %d (de 0 a %d processos).\n",
			rank, size-1);
	else
	printf( "Ola, mundo. Eu sou o processo no. %d.\n",
			rank);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
*/
