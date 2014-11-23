/*
 ====================
 JoaninhasCalorosas.c
 ====================

 ==============================================================================
 MAC0431 -- 28/11/2014 -- IME/USP, -- Prof. Marco Dimas Gubitoso
 Author      : Marcello Souza de Oliveira
 Num. USP    : 6432692
 Course      : Bacharelado em Ciencias da Computacao

 Name        : JoaninhasCalorosas.c
 Version     :
 Copyright   :
 Description : Joaninhas Calorosas: simulacao...
   ...

 Compiler         : gcc linux 4.6.3
 Compiler Options : -Wall -ansi -pedantic -O2 -U_FORTIFY_SOURCE
 Link Options     : -lm
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

 How to execute:
 $ mpicc joaninhasCalorosas.c -o joaninhasCalorosas
 $ mpirun -np 4 ./joaninhasCalorosas <params>
 ==============================================================================
*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h> /* 'linkagem' com '-lm': gcc fontes.c bibs.h -o exec -lm */
#include <time.h>

/*#include <mpi.h>*/
#include <omp.h>

/*********
 * DEFINES
 ********/

#ifndef M_PI
   #define M_PI 3.14159265358979323846  /* pi in math.h */
#endif

/**********
 * DEBUGING
 *********/

#define DEBUG    /* Habilita o modo de depuração */
#define PRINT    /* Habilita a impressao de resultados na tela */
#define CONV 9   /* Habilita a impressao da tab de convergencia */

#ifdef DEBUG
   #define MSG_DBG( fmt, ...) fprintf( stderr, "\t[MSG_DBG:] %s (%s:%d): " \
                                       fmt"\n", __FUNCTION__, __FILE__, \
                                       __LINE__, ##__VA_ARGS__)
#else
   #define MSG_DBG( fmt, ...) {;}
#endif

/*********
 * STRUCTS
 ********/

typedef enum { VAZIA, JOANINHA, CALOR, FRIO } occupation;
typedef enum { EVEN, ODD } parity;

enum hexdirection { EAST, NORTHEAST, NORTHWEST, WEST, SOUTHWEST, SOUTHEST };
typedef enum hexdirection hexdirection;

typedef struct hex {
	int energy;
	occupation elem;
} hex;
typedef hex **hexgrid;

typedef struct position {
	int row;
	int col;
} position;

typedef struct dposition {
	int row;
	int col;
	unsigned int bornscycle;
} dposition;

/*
 * possibilita testar ao final dos calculos se duas joaninhas desejam se mover
 * para o mesmo hex (devera prevalecer a com maior diferenca de energia para
 * com a energia do hex destino, conforme enunciado).
 */
typedef struct bugsmove {
	position orig;
	position dest;
} bugsmove;

struct cycledoublequeue {
    dposition *hotspos;			/* hot's source dpositions vector */
    dposition *coldspos;			/* Cold's source dpositions vector */
    int hotsbegin;
    int hotsend;
    int hotssize;
    int coldsbegin;
    int coldsend;
    int coldssize;
};
typedef struct cycledoublequeue cdqueue;

typedef struct neighbors {
	position nbs[6];
	int amount;		/* 0 - 6 vizinhos possiveis (considera-se apenas hex's VAZIA's) */
} neighbors;

/******************
 * GLOBAL VARIABLES
 *****************/

static unsigned int cycle;
static unsigned int A, L, j, s, C,
				    th_min, th_max,
					pc, nc, pf, nf,
					T, P, sem;

/************
 * PROTOTYPES
 ***********/

/***********
 * FUNCTIONS
 **********/

/* Versão da func. malloc com tratamento de erro.
 ************************************************/
void *mallocX( unsigned int nbytes) {
	void *ptr;
	ptr = malloc( nbytes);
	if (ptr == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
				 "malloc devolveu NULL! (%d bytes)\n", -1, nbytes);
		exit( EXIT_FAILURE);
	}
	return ptr;
}

/*
 * “even-r” horizontal layout.
 */
hex **createHexGrid() {
	int a, l;
	int dneighbors[2][6][2] = {
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } },
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } }
	};
	hex **H = (hex **) malloc( A * sizeof (struct hex *));

	if (H == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
			     "malloc devolveu NULL! (%d bytes)\n",
			     -1, (int) (A * sizeof (struct hex *)));
		exit( EXIT_FAILURE);
	}
	for (a = 0; a < A; ++a) {
		H[a] = (hex *) malloc( L * sizeof (struct hex));
		if (H[a] == NULL) {
			fprintf( stderr, "%d - A memoria estourou! Socorro! "
					 "malloc devolveu NULL! (%d bytes)\n",
					 -1, (int) (L * sizeof (struct hex)));
			exit( EXIT_FAILURE);
		}
		for (l = 0; l < L; ++l) {
			H[a][l].energy = 0;
			H[a][l].elem = VAZIA;
		}
	}
	return H;
}

bugsmove *createBugs( hex *H[]) {
	int i, k, r, q;
	bugsmove *p = (bugsmove *) malloc( j * sizeof (struct bugsmove));
	if (p == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
			     "malloc devolveu NULL! (%d bytes)\n",
			     -1, (int) (j * sizeof (struct position)));
		exit( EXIT_FAILURE);
	}
	for (i = 0, k = 0; i < j;) {
		r = (int) (((double) rand_r( &sem) / (RAND_MAX + 1.0)) * A);
		q = (int) (((double) rand_r( &sem) / (RAND_MAX + 1.0)) * L);
		if (H[r][q].elem == VAZIA) {
			H[r][q].elem = JOANINHA;
			p[k].dest.row = p[k].orig.row = r; p[k].dest.col = p[k].orig.col = q;
			i++; k++;
		}
	}
	return p;
}

/*
 ==============================================================================
*/

cdqueue *QUEUEinit() {
	cdqueue *cdq = malloc( sizeof (cdqueue));
	cdq->hotssize = 0.5;
	if (cdq->hotssize > pc)
		cdq->hotssize = pc;
	cdq->hotspos = (dposition *) malloc( ((int) (2 * cdq->hotssize * L * A)) * sizeof (dposition));
	if (cdq->hotspos == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
			     "malloc devolveu NULL! (%d bytes)\n",
			     -1, (int) (((int) (2 * cdq->hotssize * L * A)) * sizeof (struct dposition)));
		exit( EXIT_FAILURE);
	}
	if (cdq->coldssize > pf)
		cdq->coldssize = pf;
	cdq->coldspos = (dposition *) malloc( ((int) (2 * cdq->coldssize * L * A)) * sizeof (dposition));
	if (cdq->coldspos == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
			     "malloc devolveu NULL! (%d bytes)\n",
			     -1, (int) (((int) (2 * cdq->coldssize * L * A)) * sizeof (struct dposition)));
		exit( EXIT_FAILURE);
	}
	cdq->hotsbegin = 0;
	cdq->hotsend = 0;
	cdq->coldsbegin = 0;
	cdq->coldsend = 0;
	return cdq;
}

int QUEUEempty( cdqueue *cdq) {
   return cdq->hotsbegin == cdq->hotsend && cdq->coldsbegin == cdq->coldsend;
}

void QUEUEput( cdqueue *cdq, dposition item, occupation elem) {
	if (elem == CALOR) {
		if ((cdq->hotsend + 1 == cdq->hotsbegin || cdq->hotsend + 1 == cdq->hotssize) && cdq->hotsbegin == 0) {
			printf( "\nSocorro! Fila vai transbordar!\n");
			exit( EXIT_FAILURE);
		}
		cdq->hotspos[cdq->hotsend].row = item.row;
		cdq->hotspos[cdq->hotsend].col = item.col;
		cdq->hotspos[cdq->hotsend++].bornscycle = item.bornscycle;
		if (cdq->hotsend == cdq->hotssize)
			cdq->hotsend = 0;
	}
	else if (elem == FRIO) {
		if ((cdq->coldsend + 1 == cdq->coldsbegin || cdq->coldsend + 1 == cdq->coldssize) && cdq->coldsbegin == 0) {
			printf( "\nSocorro! Fila vai transbordar!\n");
			exit( EXIT_FAILURE);
		}
		cdq->coldspos[cdq->coldsend].row = item.row;
		cdq->coldspos[cdq->coldsend].col = item.col;
		cdq->coldspos[cdq->coldsend++].bornscycle = item.bornscycle;
		if (cdq->coldsend == cdq->coldssize)
			cdq->coldsend = 0;
	}
}

void QUEUEupdate( cdqueue *cdq) {
	while (cdq->hotspos[cdq->hotsbegin].bornscycle + nc <= cycle && cdq->hotsbegin < cdq->hotsend)
		cdq->hotsbegin++;
	while (cdq->coldspos[cdq->coldsbegin].bornscycle + nf <= cycle && cdq->coldsbegin < cdq->coldsend)
		cdq->coldsbegin++;
}

void QUEUEfree( cdqueue *cdq) {
	if (cdq->hotspos != NULL)
		free( cdq->hotspos);
	if (cdq->coldspos != NULL)
		free( cdq->coldspos);
}

/*
 ==============================================================================
*/

void generateDisturbs( hex *H[], cdqueue *cdq) {
	int a, l;
	dposition dp;
	for (a = 0; a < A; ++a)
		for (l = 0; l < L; ++l) {
			if (H[a][l].elem == VAZIA && rand_r( &sem) / (RAND_MAX + 1.0) <= pc) {
				H[a][l].elem = CALOR;
				dp.row = a; dp.col = l; dp.bornscycle = cycle;
				QUEUEput( cdq, dp, CALOR);
			}
			if (H[a][l].elem == VAZIA && rand_r( &sem) / (RAND_MAX + 1.0) <= pf) {
				H[a][l].elem = FRIO;
				dp.row = a; dp.col = l; dp.bornscycle = cycle;
				QUEUEput( cdq, dp, FRIO);
			}
		}
	return;
}

/*
 * hotspos[2*pc], coldspos[2*pf];
 */
void nextCycle( hex *H[], cdqueue *cdq) {
	cycle++;
	QUEUEupdate( cdq);
	generateDisturbs( H, cdq);
}

neighbors *neighborn( hex *H[], position p) {
	int i, k = 0,
		paritying = p.row % 2;
	neighbors *nbs = (neighbors *) malloc( sizeof (neighbors));
	int neighbors[2][6][2] = {
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } },
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } }
	};
	if (nbs == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
			     "malloc devolveu NULL! (%d bytes)\n",
			     -1, (int) (sizeof (struct neighbors)));
		exit( EXIT_FAILURE);
	}
	for (i = 0; i < 6; ++i) {
		if (H[p.row + neighbors[paritying][i][0]][p.col + neighbors[paritying][i][1]].elem == VAZIA) {
			nbs->nbs[k].row = p.row + neighbors[paritying][i][0];
			nbs->nbs[k].col = p.row + neighbors[paritying][i][1];
			nbs->amount = ++k;
		};
	}
	return nbs;
}

/*
 ==============================================================================
*/

int startSimulation( hex *H[], bugsmove *bugspos) {
	cdqueue *cdq = QUEUEinit();
	cycle = 0;

	while (cycle < T)
		nextCycle( H, cdq);

	return EXIT_SUCCESS;
}

void killBugs( bugsmove *bugspos) {
	return;
}

void destroyHexGrid( hex *H[]) {
	return;
}

/*
 ==============================================================================
*/

int main( int argc, char *argv[]) {
	bugsmove *bugspos;			/* Bug's positions vector */
	hex **H;					/* Hex Matrix (Hex Grid) */

	if (argc < 12) {
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

	scanf( " %u %u %u %u %u %u %u %u %u %u %u %u %u",
		   &L, &A, &j, &s, &C, &th_min, &th_max, &pc, &nc, &pf, &nf, &T, &P);

	H = createHexGrid( A, L);
	bugspos = createBugs( H);

	startSimulation( H, bugspos);

	killBugs( bugspos);
	destroyHexGrid( H);

	return EXIT_SUCCESS;
}



/*
 ==============================================================================
*/
	#ifdef PRINT
		/*printf( " %19.16G  %19.16G  %19.16G\n",
				t_kless[0], ab->pfunc_y( t_kless[0]), y_kless[0]);*/
	#endif


/* 	disturbpos = (position *) malloc( (A * L) * sizeof (struct position));   ***AQUI, MELHORAR ISTO: MUDAR PARA LISTA LIGADA **/

/* array[array_size - (abs(index) % array_size)] */

/*
	neighbors[EVEN][EAST]      = {  0, +1 };
	neighbors[EVEN][NORTHEAST] = { -1,  0 };
	neighbors[EVEN][NORTHWEST] = { -1, -1 };
	neighbors[EVEN][WEST]      = {  0, -1 };
	neighbors[EVEN][SOUTHWEST] = { +1, -1 };
	neighbors[EVEN][SOUTHEAST] = { +1,  0 };
	neighbors[ODD][EAST]       = {  0, +1 };
	neighbors[ODD][NORTHEAST]  = { -1, +1 };
	neighbors[ODD][NORTHWEST]  = { -1,  0 };
	neighbors[ODD][WEST]       = {  0, -1 };
	neighbors[ODD][SOUTHWEST]  = { +1,  0 };
	neighbors[ODD][SOUTHEAST]  = { +1, +1 };
*/



/*
//neighbors = [
//   [ [ 0, +1], [-1,  0], [-1, -1], [ 0, -1], [+1, -1], [+1,  0] ],
//   [ [ 0, +1], [-1, +1], [-1,  0], [ 0, -1], [+1,  0], [+1, +1] ]
//]
//parity = r & 1;
//d = neighbors[parity][direction];
//return Hex( r + d[1], q + d[0]);
*/


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
