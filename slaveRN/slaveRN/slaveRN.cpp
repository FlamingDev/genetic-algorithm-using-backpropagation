/************************************************
* Backpropagation algorithm.
*
* Training a Neural Network, or an Autoencoder.
*
* Likewise, you can increase the number of layers
* to implement a deeper structure, to follow the
* trend of "Deep Learning". But bear in mind that
* DL has lots of beautiful tricks, and merely
* making it deeper will not yield good results!!
*
* Designed by Junbo Zhao, 12/14/2013
************************************************/

#include <stdio.h>
#include <string>
//#include <afx.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#define InputN 64		// number of neurons in the input layer
#define HN 25			// number of neurons in the hidden layer
#define OutN 64			// number of neurons in the output layer
#define datanum 500		// number of training samples

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	MPI_Comm parent;
	MPI_Comm_get_parent(&parent);

	if (parent == MPI_COMM_NULL) {
		printf("Erro: Nenhum processo pai encontrado!\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Receber indivíduo e fitness
	double individuo[10];
	double fitness;
	MPI_Recv(individuo, 10, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
	MPI_Recv(&fitness, 1, MPI_DOUBLE, 0, 1, parent, MPI_STATUS_IGNORE);

	char log_msg[100];
	sprintf_s(log_msg, "Aqui eh o trabalhador %d! Recebi um individuo e fitness: %.2f\n", rank, fitness);
	
	MPI_Send(log_msg, strlen(log_msg) + 1, MPI_CHAR, 0, 2, parent);
	MPI_Finalize();
	while (1) { ; }


	double sigmoid(double);
	std::string result = "";
	char buffer[200];
	double x_out[InputN];		// input layer
	double hn_out[HN];			// hidden layer
	double y_out[OutN];         // output layer
	double y[OutN];				// expected output layer
	double w[InputN][HN];		// weights from input layer to hidden layer
	double v[HN][OutN];			// weights from hidden layer to output layer

	double deltaw[InputN][HN];
	double deltav[HN][OutN];

	double hn_delta[HN];		// delta of hidden layer
	double y_delta[OutN];		// delta of output layer
	double error;
	double errlimit = 0.001;
	double alpha = 0.1, beta = 0.1;
	int loop = 0;
	int times = 50000;
	int i, j, m;
	double max, min;
	double sumtemp;
	double errtemp;

	// training set
	struct {
		double input[InputN];
		double teach[OutN];
	}data[datanum];

	// Generate data samples
	// You can use your own data!!!
	for (m = 0; m < datanum; m++) {
		for (i = 0; i < InputN; i++)
			data[m].input[i] = (double)rand() / 32767.0;
		for (i = 0;i < OutN;i++)
			data[m].teach[i] = (double)rand() / 32767.0;
	}

	// Initializition
	for (i = 0; i < InputN; i++) {
		for (j = 0; j < HN; j++) {
			w[i][j] = ((double)rand() / 32767.0) * 2 - 1;
			deltaw[i][j] = 0;
		}
	}
	for (i = 0; i < HN; i++) {
		for (j = 0; j < OutN; j++) {
			v[i][j] = ((double)rand() / 32767.0) * 2 - 1;
			deltav[i][j] = 0;
		}
	}

	// Training
	while (loop < times) {
		loop++;
		error = 0.0;

		for (m = 0; m < datanum; m++) {
			// Feedforward
			max = 0.0;
			min = 0.0;
			for (i = 0; i < InputN; i++) {
				x_out[i] = data[m].input[i];
				if (max < x_out[i])
					max = x_out[i];
				if (min > x_out[i])
					min = x_out[i];
			}
			for (i = 0; i < InputN; i++) {
				x_out[i] = (x_out[i] - min) / (max - min);
			}

			for (i = 0; i < OutN; i++) {
				y[i] = data[m].teach[i];
			}

			for (i = 0; i < HN; i++) {
				sumtemp = 0.0;
				for (j = 0; j < InputN; j++)
					sumtemp += w[j][i] * x_out[j];
				hn_out[i] = sigmoid(sumtemp);		// sigmoid serves as the activation function
			}

			for (i = 0; i < OutN; i++) {
				sumtemp = 0.0;
				for (j = 0; j < HN; j++)
					sumtemp += v[j][i] * hn_out[j];
				y_out[i] = sigmoid(sumtemp);
			}

			// Backpropagation
			for (i = 0; i < OutN; i++) {
				errtemp = y[i] - y_out[i];
				y_delta[i] = -errtemp * sigmoid(y_out[i]) * (1.0 - sigmoid(y_out[i]));
				error += errtemp * errtemp;
			}

			for (i = 0; i < HN; i++) {
				errtemp = 0.0;
				for (j = 0; j < OutN; j++)
					errtemp += y_delta[j] * v[i][j];
				hn_delta[i] = errtemp * (1.0 + hn_out[i]) * (1.0 - hn_out[i]);
			}

			// Stochastic gradient descent
			for (i = 0; i < OutN; i++) {
				for (j = 0; j < HN; j++) {
					deltav[j][i] = alpha * deltav[j][i] + beta * y_delta[i] * hn_out[j];
					v[j][i] -= deltav[j][i];
				}
			}

			for (i = 0; i < HN; i++) {
				for (j = 0; j < InputN; j++) {
					deltaw[j][i] = alpha * deltaw[j][i] + beta * hn_delta[i] * x_out[j];
					w[j][i] -= deltaw[j][i];
				}
			}
		}

		// Global error
		error = error / 2;
		if (loop % 1000 == 0) {
			result = "Global Error = ";
			sprintf_s(buffer, "%f", error);
			result += buffer;
			result += "\r\n";
		}
		if (error < errlimit)
			break;

		printf("The %d th training, error: %f\n", loop, error);
	}

}

// sigmoid serves as avtivation function
double sigmoid(double x) {
	return(1.0 / (1.0 + exp(-x)));
}