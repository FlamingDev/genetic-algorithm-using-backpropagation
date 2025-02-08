#define GARE "GAv19"
#define NOLINHA
#define XLS
#define NOPLOT
#define NODUMP
#define CONSO
#define PASSOS 4
#define NOTESTE
#define ANSI

#ifndef ANSI

#include <condefs.h>*/
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

/* Problema: TESTAR ras ros gri sch sph  */
#define FUNCAO ras
#define MAXVAR 1000
#define SOLUCAO 0.0F

#define PI 3.14159265358979F

/* Persistencia */
#define MAXGER 250000
#define NUMCRU 100
#define TAXERR 1.0e-3
#define MAXAVA 10000000

/* Populacao*/
#define PATUAL 10
#define MAXPOP 100

/* Selecao e Cruzamento */
#define FATOR  1000

#define XALFA  0.25

/* Mutacao */
#define MNUNI  7
#define PMUTAC 5

/* Redes neurais */
#define NUM_WORKERS 4

/* ****************FUNCOES TESTE IMPLEMENTADAS *************** */

int Aval = 0;
/*************** CODIGO DAS FUNCOES TESTE ******************/


/*1.9 Ackley's Path function 10
Ackley's Path [Ack87] is a widely used multimodal test function.
function definition
f10(x)=-a·exp(-b·sqrt(1/n·sum(x(i)^2)))-exp(1/n·sum(cos(c·x(i))))+a+exp(1);
a=20; b=0.2; c=2·pi; i=1:n;
-32.768<=x(i)<=32.768.
global minimum
f(x)=0; x(i)=0, i=1:n.*/

double Ackley(double x[], int n)
{
    register int i;
    double sum;

    double a = 20.0F, b = 0.2F, c = 2.0F * PI;
    double sum1, sum2;

    for (i = 0, sum1 = 0, sum2 = 0; i < n;i++) {
        sum1 += pow(x[i], 2);
        sum2 += cos(c * x[i]);
    }
    sum = -a * exp(-b * sqrt(1.0F / n * sum1)) - exp(1.0F / n * sum2) + a + exp(1.0F);

    Aval++;

    return sum;
}

double sphere(double x[], int n)
{
    register int i;
    double sum = 100;
    for (i = 0; i < n; i++)
        sum += pow(x[i], 2);
    Aval++;

    return sum;
}
/*
 * 1st ICEO test function
 *
 * Problem 3: Shekel's Foxholes
 *
 */

double shekel(double* x, int n)
{
    double a[30][10] = {
{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
{7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
{1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
{8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
{0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
{7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
{0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
{2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
{8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
{2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
{4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
{8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
{8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
{4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
{2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
{6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
{0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
{5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
{3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
{8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
{1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
{0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
{0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
{4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
{9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
{4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500} };


    double c[] = { 0.806, 0.517, 0.1, 0.908, 0.965, 0.669, 0.524, 0.902, 0.531,
0.876, 0.462, 0.491, 0.463, 0.714, 0.352, 0.869, 0.813, 0.811, 0.828,
0.964, 0.789, 0.360, 0.369, 0.992, 0.332, 0.817, 0.632, 0.883, 0.608,
0.326 };

    register int	i, j;
    double		sp, h, result = 0.0;

    for (i = 0; i < 30; i++) {
        sp = 0.0;
        for (j = 0; j < n; j++) {
            h = x[j] - a[i][j];
            sp += h * h;
        }
        result += 1.0 / (sp + c[i]);
    }
    Aval++;

    return(-result);
}

/*
 * 1st ICEO test function
 *
 * Problem 5: The Langerman's function
 *
 */

double SqrDst(double x1[], double x2[], int n)
{
    double dist = 0.0, d;
    int i;

    dist = 0;

    for (i = 0; i < n; i++)
    {
        d = x1[i] - x2[i];
        dist += d * d;
    }

    return (dist);
}


double langerman(double 	x[], int		n)
{
    double a[30][10] = {
{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
{7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
{1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
{8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
{0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
{7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
{0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
{2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
{8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
{2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
{4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
{8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
{8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
{4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
{2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
{6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
{0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
{5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
{3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
{8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
{1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
{0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
{0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
{4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
{9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
{4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500} };


    double c[] = {
    0.806,
0.517,
0.1,
0.908,
0.965,
0.669,
0.524,
0.902,
0.531,
0.876,
0.462,
0.491,
0.463,
0.714,
0.352,
0.869,
0.813,
0.811,
0.828,
0.964,
0.789,
0.360,
0.369,
0.992,
0.332,
0.817,
0.632,
0.883,
0.608,
0.326 };
    //        double min = 500;

    register int 	i;
    double 		Sum, exp(double),
        dist;
    Sum = 0;
    for (i = 0; i < 5; i++)
    {
        dist = SqrDst(x, a[i], n);
        Sum -= c[i] * (exp(-dist / PI) * cos(PI * dist));
    }
    Aval++;

    return (Sum);
}


/*      ICEO Michalewitz */

double Micha2(double x[], int n)
{
    double m = 10.0F;
    double u = 0;
    register int i;

    for (i = 0; i < n; i++) {
        u = u - sin(x[i]) * pow(sin(i * pow(x[i], 2) / PI), 2 * m);
    }

    return u;


}
double Micha(double y[], int n)
{
    double x[MAXVAR];
    double   u, temp, m = 10.0F;
    double cost, sint;
    int     i;

    u = 0;

    cost = cos(PI / 6);
    sint = sin(PI / 6);

    for (i = 0; i < n - 1; i += 2)
    {
        x[i] = y[i] * cost - y[i + 1] * sint;
        x[i + 1] = y[i] * sint + y[i + 1] * cost;
    }

    if (i == n - 1)
        x[i] = y[i];

    for (i = 0;i < n;i++)
        u = u + sin(x[i])
        * pow(sin((i + 1) * x[i] * x[i] / PI), 2.0 * m);

    Aval++;

    return(-u);
}

/*      function definition
        f12(x)=-sum(sin(x(i))·(sin(i·x(i)^2/pi))^(2·m)), i=1:n, m=10;
        0<=x(i)<=pi.
        global minimum
        f(x)=-4.687 (n=5); x(i)=???, i=1:n.
        f(x)=-9.66 (n=10); x(i)=???, i=1:n.
*/
double fun12(double x[], int n)
{
    register int i;
    int m = 10;
    double sum = 0;

    for (i = 0; i < n; i++)
        sum += sin(x[i]) * pow(sin((i + 1) * pow(x[i], 2) / PI), (2 * m));
    Aval++;

    return -sum;
}

/*      function definition
        fEaso(x1,x2)=-cos(x1)·cos(x2)·exp(-((x1-pi)^2+(x2-pi)^2));
        -100<=x(i)<=100, i=1:2.
        global minimum
        f(x1,x2)=-1; (x1,x2)=(pi,pi).
*/
double fEaso(double x[], int n)
{
    register int i;
    double sum;

    sum = -cos(x[1]) * cos(x[2]) * exp(-(pow((x[1] - PI), 2) + pow((x[2] - PI), 2)));

    Aval++;

    return sum;
}
/*      function definition
        fGold(x1,x2)=[1+(x1+x2+1)^2·(19-14·x1+3·x1^2-14·x2+6·x1·x2+3·x2^2)]·[30+(2·x1-3·x2)^2·(18-32·x1+12·x1^2+48·x2-36·x1·x2+27·x2^2)];
        -2<=x(i)<=2, i=1:2.
        global minimum
        f(x[1],x[2])=3; (x[1],x[2])=(0,-1).
*/
double fGold(double x[], int n)
{
    register int i;
    double sum;

    sum = (1 + pow((x[1] + x[2] + 1), 2) * (19 - 14 * x[1] + 3 * pow(x[1], 2) - 14 * x[2] + 6 * x[1] * x[2] + 3 * pow(x[2], 2))) * (30 + pow((2 * x[1] - 3 * x[2]), 2) * (18 - 32 * x[1] + 12 * pow(x[1], 2) + 48 * x[2] - 36 * x[1] * x[2] + 27 * pow(x[2], 2)));

    Aval++;

    return sum;
}

double schwefel7(double x[], int n)
{
    register int i;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += -x[i] * sin(sqrt(fabs(x[i])));
    }
    Aval++;

    return sum;
}

double rastringin(double x[], int n)
{
    register int i;
    double sum;
    for (i = 0, sum = 3 * n; i < n; i++)
    {
        sum = sum + pow(x[i], 2) - 3 * cos(2 * PI * x[i]);
    }
    Aval++;

    return sum;
}
double griewangk(double x[], int n)
{
    register int i;
    double sum, pro;
    for (i = 1, sum = 0, pro = 1; i <= n; i++) {
        sum = sum + pow(x[i - 1], 2) / 4000;
        pro = pro * cos(x[i - 1] / sqrt(i));
    }
    Aval++;


    return (1 + sum + pro);
}

double rosenbrock(double x[], int n)
{
    register int i;
    double sum = 0;
    for (i = 0; i < n - 1; i++)
    {
        sum += 100 * pow(x[i + 1] - pow(x[i], 2), 2) + pow(1 - x[i], 2);
    }
    Aval++;

    return sum;
}


/* *************** ESTRUTURA DAS FUNCOES TESTE ***************** */

enum { f12 = 0, gol, fea, sph, ros, sch, ras, gri, mic, lan, she, ack } enumFunc;
static double (*funccod[])(double[], int n) =
{ fun12, fGold, fEaso, sphere, rosenbrock, schwefel7, rastringin, griewangk, Micha, langerman, shekel, Ackley };
struct DadosFuncoes {
    char    nom[5];
    double  inf;
    double  sup;
}limFixos[] = { {"F12", 0.0F, PI},{"Gol",-2.0F,2.0F},{"Fea",-100.0F,100.0F},
                {"Sph",-5.12F, 5.12F},{"Ros",-2.048F, 2.048F},{"Sch",-500.0F, 500.0F},
                {"Ras",-5.12,5.12}, {"Gri",-600.0, 600.0}, {"Mic",-0.0F, PI},
                {"Lan",0.0F, 10.0F}, {"She",-1000.0F, 1000.0F}, {"ack",-32.768F, 32.768F} };

int funcao = FUNCAO;

typedef struct {
    double          var[MAXVAR];
    double          fit;
}Cromossomo;

typedef struct {
    Cromossomo      indiv[MAXPOP];
    double          sumFit;
    int             tamPop;
    int             tamInd;
    int             melhor;
    int             pior;
    int             numMuta;
    int             iguais;
    int             gerMelhor;
} Populacao;



float randgen(float fLlim, float fUlim)
{
    float fRandomVal;

    fRandomVal = rand() % 101 / 100.;        // rand entre 0 e 1

    return(fLlim + (float)(fRandomVal * (fUlim - fLlim)));

}

void nu_mutate(double* indiv, int tamind, int tampop, int ger, int expo)
{
    int i, j;
    float fRandVal, fFactor;
    float fNewt, fNewT;
    int   iExponent, iIndex;

    fRandVal = (rand() % 101 / 100.);
    iIndex = rand() % tamind;
    /* pick either the max or min. limit */
    if (fRandVal < 0.5)		/* lower */
    {
        fNewt = ger;
        fNewT = MAXGER;
        //                iExponent = (int) (6*randgen(0.0, 1.0)) + 1;
        fRandVal = (rand() % 101 / 100.);
        fFactor = pow((1.0F - (fNewt / fNewT)), expo) * fRandVal;

        if (fFactor < TAXERR) fFactor = TAXERR;
        fFactor = fFactor * (indiv[iIndex] - limFixos[funcao].inf);
        indiv[iIndex] = indiv[iIndex] - fFactor;
    }
    else
    {
        fNewt = ger;
        fNewT = MAXGER;
        //                iExponent = (int) (6*randgen(0.0, 1.0)) + 1;
        fRandVal = (rand() % 101 / 100.);

        fFactor = pow((1.0F - (fNewt / fNewT)), expo) * fRandVal;

        if (fFactor < TAXERR) fFactor = TAXERR;
        fFactor = fFactor * (limFixos[funcao].sup - indiv[iIndex]);
        indiv[iIndex] = indiv[iIndex] + fFactor;
    }
}


void IniciaPop(Populacao* p, int m, int n, int nfun) {
    int i, j, pior, melhor;
    double soma, fit;

    for (i = 0, soma = 0, pior = 0, melhor = 0; i < m; i++)
    {
        for (j = 0; j < n; j++) {
            p->indiv[i].var[j] = (double)randgen(limFixos[nfun].inf, limFixos[nfun].sup);
        }
        fit = funccod[nfun](p->indiv[i].var, n);
        p->indiv[i].fit = fit;
        if (fit > p->indiv[pior].fit) pior = i;
        if (fit < p->indiv[melhor].fit) melhor = i;
        soma += (fit); //fabs
    }
    p->tamPop = m;
    p->tamInd = n;
    p->sumFit = soma;
    p->melhor = melhor;
    p->pior = pior;
    p->numMuta = 0;
    p->iguais = 0;
}
void IniciaPSc(Populacao* p, int m, int n, int nfun) {
    int i, j, k;
    int pior, melhor, ipp;
    float inf, sup, pas, par, npa;
    double soma, fit;

    npa = log(m) / log(2);
    pas = (limFixos[nfun].sup - limFixos[nfun].inf) / npa;
    inf = limFixos[nfun].inf;
    sup = ((inf + pas) < limFixos[nfun].sup ? (inf + pas) : limFixos[nfun].sup);
    ipp = m / npa;

    for (i = 0, soma = 0, pior = 0, melhor = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            p->indiv[i].var[j] = (double)randgen(inf, sup);
        fit = funccod[nfun](p->indiv[i].var, n);
        p->indiv[i].fit = fit;
        if (fit > p->indiv[pior].fit) pior = i;
        if (fit < p->indiv[melhor].fit) melhor = i;
        soma += (fit); //fabs
        if ((i + 1) % ipp == 0) {
            inf = sup;
            sup = ((inf + pas) < limFixos[nfun].sup ? (inf + pas) : limFixos[nfun].sup);
        }
    }
    p->tamPop = m;
    p->tamInd = n;
    p->sumFit = soma;
    p->melhor = melhor;
    p->pior = pior;
    p->numMuta = 0;
    p->iguais = 0;
}

int  Selecao(Populacao* p)
{
    int   i;
    double val, spin_val;

    val = 0.0;
    spin_val = (rand() % 101) / 100. * p->sumFit;
    i = rand() % p->tamPop;

    do {
        i = (i < p->tamPop - 1) ? i + 1 : 0;
        val = val + fabs(FATOR * p->sumFit) / (1 + p->indiv[i].fit - p->indiv[p->melhor].fit);
    } while (val < spin_val);
    return i;        // posição do estouro
}

void CruzaBlend(Populacao* p, int pai, int mae, int filho, float alfa)
{
    float a, b, r;
    int i;
    double ajuste, nureal, interv;
    int nuinte, voltas;


    a = -alfa;
    b = 1 + alfa;
    r = rand() % 101 / 100.;
    r = a + r * (b - a);

    for (i = 0;i < p->tamInd; i++) {
        p->indiv[filho].var[i] =
            p->indiv[pai].var[i] + r * (p->indiv[mae].var[i] - p->indiv[pai].var[i]);
        if (p->indiv[filho].var[i] < limFixos[FUNCAO].inf) {
            interv = limFixos[FUNCAO].sup - limFixos[FUNCAO].inf;
            nureal = (limFixos[FUNCAO].inf - p->indiv[filho].var[i]) / interv;
            nuinte = nureal;
            ajuste = (nureal - nuinte) * interv;
            voltas = nuinte % 2;
            p->indiv[filho].var[i] = limFixos[FUNCAO].inf + voltas * interv - (2 * voltas - 1) * ajuste;
        }
        else if (p->indiv[filho].var[i] > limFixos[FUNCAO].sup) {
            interv = limFixos[FUNCAO].sup - limFixos[FUNCAO].inf;
            nureal = (p->indiv[filho].var[i] - limFixos[FUNCAO].sup) / interv;
            nuinte = nureal;
            ajuste = (nureal - nuinte) * interv;
            voltas = 1 - nuinte % 2;
            p->indiv[filho].var[i] = limFixos[FUNCAO].inf + voltas * interv - (2 * voltas - 1) * ajuste;
        }
    }

}


void AtualizaPop(Populacao* p, int pos, double fit, int ger)
{
    int i, j, max;

    p->sumFit -= (p->indiv[pos].fit); //fabs
    p->sumFit += (fit);               // fabs
    p->indiv[pos].fit = fit;

    if (fit < p->indiv[p->melhor].fit) {
        p->melhor = pos;
        p->gerMelhor = ger;
    }
    /* ***** procura um outro pior ******** */
    max = PATUAL * p->tamPop / 100;
    //        p->pior=p->melhor == p->tamPop ? p->melhor - 1: p->melhor +1;
    for (i = 0; i < max; i++) {
        j = rand() % p->tamPop;
        if (j == pos || j == p->melhor)
            continue;
        if (p->indiv[j].fit > p->indiv[p->pior].fit) {
            p->pior = j;
            i += PATUAL;
        }
    }

}

/****/



/*******************************************************************/

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm intercomm;
    const char* path = "C:/Users/m1r9b/source/repos/genetic-algorithm-using-backpropagation/slaveRN/x64/Debug/slaveRN.exe";
    int workers_ids[NUM_WORKERS];

    MPI_Comm_spawn(path, MPI_ARGV_NULL,NUM_WORKERS, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &intercomm, workers_ids);

    // individuo e fitness teste
    double individuo[10] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    double fitness = 0.5;


    for (int i = 0; i < NUM_WORKERS; i++) { // envia individuo e fitness para todas RNs spawnadas
        MPI_Send(individuo, 10, MPI_DOUBLE, i, 0, intercomm); // tag 0 para individuo
        MPI_Send(&fitness, 1, MPI_DOUBLE, i, 1, intercomm); // tag 1 para fitness
    }

    printf("Mestre enviou os dados para os trabalhadores.\n");

    char msg[100];

    for (int i = 0; i < NUM_WORKERS; i++) {
        MPI_Recv(msg, 100, MPI_CHAR, i, 2, intercomm, MPI_STATUS_IGNORE); // tag 2 para mensagens
        printf("%s", msg);
        fflush(stdout);
    }

    clock_t start, end;

    FILE* saida = NULL, * arq = NULL;
    char nomarq[50];

    int MaxIt = PASSOS, numGeracoes = MAXGER, numCruza;
    float pMuta = PMUTAC, xalfa = XALFA;
    double erro;

    Populacao P;
    int pa1, pa2;
    double fit, dvp, med;

    int i, j;

    int semente;


#ifdef LINHA
    funcao = atoi(argv[2]);
    strcpy(nomarq, GARE);
    strcat(nomarq, limFixos[funcao].nom);
    strcat(&nomarq[strlen(limFixos[funcao].nom)], argv[1]);
    if (!(saida = fopen(nomarq, "w"))) {
        perror("");
        exit(-1);
    }

    funcao = atoi(argv[2]);
    semente = atoi(argv[3]);
#endif

//#ifdef NOLINHA
    saida = stdout;
//#endif

    // randomico ou não 
    srand((unsigned)time(0));

    IniciaPop(&P, MAXPOP, MAXVAR, funcao);
    erro = (double)P.indiv[P.melhor].fit - SOLUCAO;

    printf("\n***   (%s) por ACMO/CAP/LAC/INPE    ***", GARE);
    printf("\n      Evoluindo %s com %d vars", limFixos[funcao].nom, MAXVAR);
    start = clock();
    while (numGeracoes-- && Aval < MAXAVA && erro >(double) 0.0F) {
        numCruza = NUMCRU;
        while (numCruza--) {
            pa1 = Selecao(&P);
            pa2 = Selecao(&P);
            if (pa1 != pa2)
                CruzaBlend(&P, pa1, pa2, P.pior, xalfa);
            else
                P.iguais++;
            if (pMuta > rand() % 100) {
                nu_mutate(P.indiv[P.pior].var, P.tamInd, P.tamPop, MAXGER - numGeracoes, MNUNI);
                P.numMuta++;
            }
            fit = funccod[funcao](P.indiv[P.pior].var, P.tamInd);
            AtualizaPop(&P, P.pior, fit, MAXGER - numGeracoes);
        }
        erro = (double)P.indiv[P.melhor].fit - SOLUCAO;
    }

    end = clock();

    // *******************RESULTADOS ********************

    med = P.sumFit / P.tamPop;
    dvp = 0;
    for (i = 0;i < P.tamPop;i++)
        dvp += pow(P.indiv[i].fit - med, 2);
    dvp = sqrt(dvp / (P.tamPop - 1));

#ifdef XLS

    fprintf(saida, "\n%s)Min = %.10f; Aval= %d; Tempo = %.4f; Med = %.4f; Dpd = %.4f; Ger = %d",
        nomarq, P.indiv[P.melhor].fit, Aval, (double)(end - start) / 118, med, dvp, P.gerMelhor);

#endif

#ifdef CONSO
    printf("\n\t Variaveis ...");
    for (i = 0;i < P.tamInd;i++) {
        printf("\n\t\t%.6f", P.indiv[P.melhor].var[i]);
    }

    printf("\n\t Minimo = %.16f; \n\t Na ger = %d; \n\t mutacoes = %d \n\t ; aval = %d\n\t; (%.4f, %.4f)\n",
        P.indiv[P.melhor].fit, P.gerMelhor, P.numMuta, Aval, med, dvp);
    printf("\n\t em tempo = %.4f, ", float((end - start) / CLOCKS_PER_SEC));
    fflush(stdout);
#endif

#ifdef DUMP

    fprintf(saida, "\nD U M P **************");
    fprintf(saida, "\n\tPopulacao{Tpop=%d, Tind=%d, Som=%.4f, X*=%d, X =%d, #m=%d, #i=%d, (%.4f, %.4f)",
        P.tamPop, P.tamInd, P.sumFit, P.melhor, P.pior, P.numMuta, P.iguais, med, dvp);

    for (i = 0;i < P.tamPop;i++) {
        fprintf(saida, "\n\tIndividuo (%d) = %.10f", i, P.indiv[i].fit);
        for (j = 0;j < P.tamInd; j++)
            fprintf(saida, "\n\t\t%.4f", P.indiv[i].var[j]);
    }
    fprintf(saida, "\n\tMelhor=%.10f,Media=%.10f, Desvio=%.10f", P.indiv[P.melhor].fit, med, dvp);
#endif
    MPI_Finalize();
    return 0;
}