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

#define DEBUG   	/* Habilita o modo de depuração */
#define PRINT   	/* Habilita a impressao de resultados na tela */
#define CONV 9  	/* Habilita a impressao da tab de convergencia */
#define NOCALC -1

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

typedef enum { VACANT = 0, LADYBUG, HOT, COLD } occupation;
typedef enum { D_HOT = 0, D_COLD = 1 } disType;
typedef enum { EVEN, ODD } parity;
/*
typedef occupation hex;
typedef occupation **Hex;
*/
typedef struct hex {
	int index;			/* indice para o vetor de joaninhas */
	occupation elem;
} hex;
typedef hex **hexgrid;

typedef struct position {
	int row;
	int col;
} position;

typedef struct dposition {
	position pos;
	float energy;
	unsigned int bornscycle;
} dposition;

typedef struct eposition {
	position pos;
	float energy;
} eposition;

/*
 * possibilita testar ao final dos calculos se duas joaninhas desejam se mover
 * para o mesmo hex (devera prevalecer a com maior diferenca de energia para
 * com a energia do hex destino, conforme enunciado).
 * amount: 0 - 6 vizinhos possiveis (considera-se apenas hex's VAZIA's).
 */
typedef struct bug {
	eposition orig;		/* ATENCAO AQUI, posso ter de mallocar orig e dest */
	eposition dest;
	eposition neighbor[6];
	int amount;
	int movesflag;
} bug;

typedef struct entData {
	unsigned int A, L, j, s, C, th_min, th_max, nc, nf, T, P, sem;
	float pc, pf;	/* prob em decimal */
} entData;

/*
 * hotspos: hot's source dpositions array (aceita ate (pc+TOL)*A*L hot's
 * sources e expansivel mais (A*L-j)-(pc+TOL)*A*L, totalizando A*L-j).
 * coldspos: Cold's source dpositions array (aceita ate (pf+TOL)*A*L cold's
 * sources e expansivel mais (A*L-j)-(pf+TOL)*A*L, totalizando A*L-j).
 */
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

/*static unsigned int cycle;*/

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
hex **createHexGrid( const entData *data) {
	int a, l;
	hex **H = (hex **) malloc( data->A * sizeof (hex *));
	if (H == NULL) {
		fprintf( stderr, "%d - A memoria estourou! Socorro! "
				 "malloc devolveu NULL! (%d bytes)\n",
				 -1, (int) (data->A * sizeof (hex *)));
		exit( EXIT_FAILURE);
	}
	for (a = 0; a < data->A; ++a) {
		H[a] = (hex *) malloc( data->L * sizeof (hex));
		if (H[a] == NULL) {
			fprintf( stderr, "%d - A memoria estourou! Socorro! "
					 "malloc devolveu NULL! (%d bytes)\n",
					 -1, (int) (data->L * sizeof (hex)));
			exit( EXIT_FAILURE);
		}
		for (l = 0; l < data->L; ++l) {
			H[a][l].elem = VACANT;
			H[a][l].index = data->j;
		}
	}
	return H;
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
		r = (int) (((double) rand_r( &data->sem) / (RAND_MAX + 1.0)) * data->A);
		q = (int) (((double) rand_r( &data->sem) / (RAND_MAX + 1.0)) * data->L);
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

void showBugs( bug ladybug[], int nbugs) {
	int b;
	for (b = 0; b < nbugs; ++b) {
		printf( "H( %d, %d)\n", ladybug[b].orig.pos.row, ladybug[b].orig.pos.col);
MSG_DBG( "H( %d, %d)", ladybug[b].orig.pos.row, ladybug[b].orig.pos.col);

	}
	printf( "\n");
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
	for (I = D_HOT; I <= D_COLD; ++I) {
		dis[I].size = ceil( min( dis[I].prob + TOL, 1.0) * data->L * data->A);
		dis[I].begin = 0;
		dis[I].end = 0;
		dis[I].pos = (dposition *) malloc( dis[I].size * sizeof (dposition));
		if (dis[I].pos == NULL) {
			fprintf( stderr, "%d - A memoria estourou! Socorro! "
					 "malloc devolveu NULL! (%d bytes)\n",
					 -1, (int) (dis[I].size * sizeof (struct dposition)));
			exit( EXIT_FAILURE);
		}
	}
}

/*
 *
 */
void QUEUEput( disturbs dis[], disType elem, int r_row, int q_col, unsigned int bornscycle) {
	/* verifica se dis[elem] esta cheio e se estiver, aloca novas posicoes */
	if (dis[elem].end + 1 == dis[elem].begin || (dis[elem].end + 1 == dis[elem].size && dis[elem].begin == 0)) {
		printf( "\nSocorro! Fila vai transbordar!\n");
		exit( EXIT_FAILURE);
	}
	dis[elem].pos[dis[elem].end].pos.row = r_row;
	dis[elem].pos[dis[elem].end].pos.col = q_col;
	dis[elem].pos[dis[elem].end].energy = 0.0;
	dis[elem].pos[dis[elem].end].bornscycle = bornscycle;
	(dis[elem].end + 1 == dis[elem].size) ? dis[elem].end = 0 : dis[elem].end++;

if (elem == D_HOT)
MSG_DBG( "HOT (%d, %d) : begin = %d : end = %d : size = %d", r_row, q_col, dis[elem].begin, dis[elem].end, dis[elem].size);
else
MSG_DBG( "COLD (%d, %d) : begin = %d : end = %d : size = %d", r_row, q_col, dis[elem].begin, dis[elem].end, dis[elem].size);

}

void QUEUEupdate( hex *H[], disturbs dis[], unsigned int cycle) {
	int I;
	for (I = D_HOT; I < D_COLD; ++I) {
		if (dis[I].begin != dis[I].end)
			while (dis[I].pos[dis[I].begin].bornscycle + dis[I].ncycle <= cycle && dis[I].begin != dis[I].end) {
MSG_DBG( "bornscycle = %u : cycle = %u", dis[I].pos[dis[I].begin].bornscycle, cycle);

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
	int a, l;
	for (a = 0; a < data->A; ++a)
		for (l = 0; l < data->L; ++l) {
			if (H[a][l].elem == VACANT) {
				if (rand_r( &data->sem) / (RAND_MAX + 1.0) <= data->pc) {
					H[a][l].elem = HOT;
					QUEUEput( dis, D_HOT, a, l, cycle);
				}
				else if (rand_r( &data->sem) / (RAND_MAX + 1.0) <= data->pf) {
					H[a][l].elem = COLD;
					QUEUEput( dis, D_COLD, a, l, cycle);
				}
			}
		}
	return;
}

/*
neighbors *neighborn( hex *H[], int L, int A, position p) {
*/
void neighbornsBugs( hex *H[], int L, int A, bug ladybug[], int nbugs) {
	int i, j, k,
		paritying,
		row, col;
	int directions[2][6][2] = {
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } },
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } }
	};
	for (j = 0; j < nbugs; ++j) {
		paritying = ladybug[j].orig.pos.row & 1;
		for (i = 0, k = 0; i < 6; ++i) {
			row = ladybug[j].orig.pos.row + directions[paritying][i][0];
			col = ladybug[j].orig.pos.col + directions[paritying][i][1];
			if (row > 0 && row < A && col > 0 && col < L && H[row][col].elem == VACANT) {
				ladybug[j].neighbor[k].pos.row = row;
				ladybug[j].neighbor[k].pos.col = col;
				ladybug[j].neighbor[k].energy = 0.0;
				ladybug[j].amount = ++k;
			}
		}
	}
}

void bugsMove( hex *H[], const entData *data, bug ladybug[]) {
	int i, b1, b2,
		paritying,
		row, col;
	int directions[2][6][2] = {
		{ {  0, +1 }, { -1,  0 }, { -1, -1 }, {  0, -1 }, { +1, -1 }, { +1,  0 } },
		{ {  0, +1 }, { -1, +1 }, { -1,  0 }, {  0, -1 }, { +1,  0 }, { +1, +1 } }
	};
	for (b1 = 0; b1 < data->j; ++b1) {
		if (ladybug[b1].movesflag == 1) {
			paritying = ladybug[b1].dest.pos.row & 1;
			for (i = 0, b2 = b1; i < 6; ++i) {
				row = ladybug[b2].dest.pos.row + directions[paritying][i][0];
				col = ladybug[b2].dest.pos.col + directions[paritying][i][1];
				if (ladybug[H[row][col].index].movesflag == 1 && row > 0 &&
				  row < data->A && col > 0 && col < data->L && H[row][col].elem == LADYBUG) {
					if (ladybug[H[row][col].index].dest.pos.row == ladybug[b2].dest.pos.row &&
					  ladybug[H[row][col].index].dest.pos.col == ladybug[b2].dest.pos.col) {
						if (abs( ladybug[H[row][col].index].orig.energy - ladybug[b2].dest.energy)
						  > abs( ladybug[b2].orig.energy - ladybug[b2].dest.energy)) {
							ladybug[b2].movesflag = 0;
							b2 = H[row][col].index;
						}
						else if (abs( ladybug[H[row][col].index].orig.energy - ladybug[b2].dest.energy)
						  < abs( ladybug[b2].orig.energy - ladybug[b2].dest.energy)) {
							ladybug[H[row][col].index].movesflag = 0;
						}
						else {
							ladybug[b2].movesflag = 0;
							ladybug[H[row][col].index].movesflag = 0;
							if (++i < 6) {
								do {
									row = ladybug[b2].dest.pos.row + directions[paritying][i][0];
									col = ladybug[b2].dest.pos.col + directions[paritying][i][1];
									b2 = H[row][col].index;
								} while (++i < 6 && ladybug[b2].movesflag == 0);
							}
						}
					}
				}
			}
		}
	}
	for (b1 = 0; b1 < data->j; ++b1) {
		if (ladybug[b1].movesflag == 1) {
			H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].elem = VACANT;
			H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].index = NOCALC;
			ladybug[b1].orig.pos.row = ladybug[b1].dest.pos.row;
			ladybug[b1].orig.pos.col = ladybug[b1].dest.pos.col;
			H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].elem = LADYBUG;
			H[ladybug[b1].orig.pos.row][ladybug[b1].orig.pos.col].index = b1;
			ladybug[b1].movesflag = 0;
		}
	}
	showBugs( ladybug, data->j);
}

/*
 * Apos o somatorio das partialEnergy, deve-se multiplicar pela cte C.
 * Quando uma joaninha se move, deve-se atualizar sua energia para 0.
 */
float partialEnergy( position p, position q, occupation q_elem) {
	float x, y;
	if (q_elem == VACANT)
		return 0.0;
	else {
		x = abs( p.row - q.row);
		y = abs( p.col - q.col) - (float) ((p.row + q.row) & 1) * 0.50;
	}
	/* JOANINHA ou CALOR */
	if (q_elem != COLD)
		return 1 / ( 0.75 * x * x + y * y);
	else
		return -1 / ( 0.75 * x * x + y * y);
}

/*
 * se energy != 0.0, suponho que seja um hex cuja energia ja foi calculada, e
 * apenas a devolvo, sem recalcular. C.c., mesmo que ja tenha sido calculada
 * (fontes que se compensem), calculo.
 * Importante: sempre que mover uma joaninha, deve-se zerar a energia de sua
 * posicao, de todos os seus vizinhos e ainda o amount (numero de vizinhos).
 */
float disturbsEnergy( eposition *orig, disturbs dis[], unsigned int C) {
	int c, I;
	if (orig->energy == 0.0) {
		/*#pragma omp parallel for private (x) reduction (+: sum)*/
		for (I = D_HOT; I < D_COLD; ++I) {
			/*#pragma omp parallel for private (x) reduction (+: sum)*/
			for (c = dis[I].begin; c != dis[I].end;) {
				orig->energy += partialEnergy( orig->pos, dis[I].pos[c].pos, HOT);
				(c + 1 == dis[I].size) ? c = 0 : c++;
			}
		}
	}
MSG_DBG( "H( %d, %d) : energy = %f", orig->pos.row, orig->pos.col, orig->energy);

	return orig->energy;
}

void bugsEnergy( hex *H[], disturbs dis[], bug ladybug[], const entData *data) {
	int b1, b2, k, minIndex, maxIndex;
	float minE, maxE, tmpE;
	/*#pragma omp parallel for private (x) reduction (+: sum)*/
	for (b1 = 0; b1 < data->j; ++b1) {
		neighbornsBugs( H, data->L, data->A, ladybug, data->j);
		/* if (ladybug[b1].orig.energy == 0.0) */
		(void) disturbsEnergy( &ladybug[b1].orig, dis, data->C);
		/*#pragma omp parallel*/
		/*#pragma omp parallel for private (x) reduction (+: sum)*/
		for (b2 = 0; b2 < data->j; ++b2) {
			if (b2 != b1) {
				ladybug[b1].orig.energy += partialEnergy( ladybug[b1].orig.pos, ladybug[b2].orig.pos, LADYBUG);
			}
		}
		ladybug[b1].orig.energy *= data->C;
MSG_DBG( "H( %d, %d) : energy = %f", ladybug[b1].orig.pos.row, ladybug[b1].orig.pos.col, ladybug[b1].orig.energy);

		if (ladybug[b1].amount > 0) {
			minE = maxE = disturbsEnergy( &ladybug[b1].neighbor[0], dis, data->C);
			minIndex = maxIndex = 0;
MSG_DBG( "H( %d, %d) : energy = %f", ladybug[b1].neighbor[0].pos.row, ladybug[b1].neighbor[0].pos.col, ladybug[b1].neighbor[0].energy);

			for (k = 1; k < ladybug[b1].amount; ++k) {
				tmpE = disturbsEnergy( &ladybug[b1].neighbor[k], dis, data->C);
MSG_DBG( "H( %d, %d) : energy = %f", ladybug[b1].neighbor[k].pos.row, ladybug[b1].neighbor[k].pos.col, ladybug[b1].neighbor[k].energy);

				if (tmpE < minE) {
					minE = tmpE;
					minIndex = k;
				}
				if (tmpE > maxE) {
					maxE = tmpE;
					minIndex = k;
				}
			}
			if (ladybug[b1].orig.energy < data->th_min) { /* procura fogo */
				ladybug[b1].dest.pos.row = ladybug[b1].neighbor[maxIndex].pos.row;
				ladybug[b1].dest.pos.col = ladybug[b1].neighbor[maxIndex].pos.col;
				ladybug[b1].movesflag = 1;
			}
			else if (ladybug[b1].orig.energy > data->th_min) { /* procura gelo */
				ladybug[b1].dest.pos.row = ladybug[b1].neighbor[minIndex].pos.row;
				ladybug[b1].dest.pos.col = ladybug[b1].neighbor[minIndex].pos.col;
				ladybug[b1].movesflag = 1;
			}
		}
	}
	bugsMove( H, data, ladybug);
}

void nextCycle( hex *H[], disturbs dis[], entData *data, unsigned int cycle) {
	QUEUEupdate( H, dis, cycle);
	generateDisturbs( H, dis, data, cycle);
}

/*
 ==============================================================================
*/

int startSimulation( hex *H[], bug ladybug[], entData *data) {
	int cycle = 0;
	disturbs dis[2];

	QUEUEinit( dis, data);

	for (cycle = 1; cycle <= data->T; ++cycle) {
		nextCycle( H, dis, data, cycle);
		bugsEnergy( H, dis, ladybug, data);
	}

	return EXIT_SUCCESS;
}

void killBugs( bug *ladybug) {
	return;
}

void destroyHexGrid( hex *H[]) {
	return;
}

/*
 ==============================================================================
*/

int main( int argc, char *argv[]) {
	int i;
	entData data;
	bug *ladybug;			/* Bug's positions vector */
	hex **H;				/* Hex Matrix (Hex Grid) */
	position p, q;
/*
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

	scanf( " %u %u %u %u %u %u %u %f %u %f %u %u %u",
		   &L, &A, &j, &s, &C, &th_min, &th_max, &pc, &nc, &pf, &nf, &T, &P);
*/
	/* verificar: th_min < th_max */
	data.L = 10;
	data.A = 10;
	data.j = 3;
	data.s = 0;
	data.C = 1;
	data.th_min = 2.0;
	data.th_max = 2.5;
	data.pc = 0.3;
	data.nc = 2;
	data.pf = 0.3;
	data.nf = 1;
	data.T = 10;
	data.P = 1;

	H = createHexGrid( &data);
	ladybug = createBugs( H, &data);
/*
for (i = 0; i < data.j; ++i)
MSG_DBG( "H( %d, %d) : energy = %f", ladybug[i].orig.pos.row, ladybug[i].orig.pos.col, ladybug[i].orig.energy);
*/
	showBugs( ladybug, data.j);
	startSimulation( H, ladybug, &data);
	showBugs( ladybug, data.j);
/*
	killBugs( ladybug);
	destroyHexGrid( H);
*/
	p.row = 1;
	p.col = 1;
	q.row = 3;
	q.col = 1;
/*
	printf( "\nD(%d,%d)->(%d,%d)= %f\n", p.row, p.col, q.row, q.col, hexDist( p, q));
	printf( "\n%d\n", (6 & 1));
	printf( "\n%d\n", (3 / 2));
*/

	return EXIT_SUCCESS;
}



/*
 ==============================================================================
*/
	#ifdef PRINT
		/*printf( " %19.16G  %19.16G  %19.16G\n",
				t_kless[0], ab->pfunc_y( t_kless[0]), y_kless[0]);*/
	#endif

/*MSG_DBG( "t_0=%f\ty_0=%f\th=%f\tt_n=%f\n", ab->t_0, ab->y_0, h, ab->t_n);*/

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
