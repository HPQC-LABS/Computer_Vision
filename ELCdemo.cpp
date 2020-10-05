/*//////////////////////////////////////////////////////////////////////////////////////////////////
///  ELCdemo.cpp    Clique Reduction by Excludable Local Configuration: Denoising Demo
///  Version 1.04         September 12th, 2014
////////////////////////////////////////////////////////////////////////////////////////////////////

Copyright 2014 Hiroshi Ishikawa. All rights reserved.
This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Hiroshi Ishikawa).

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For the latest version, check: http://www.f.waseda.jp/hfs/indexE.html

////////////////////////////////////////////////////////////////////////////////////////////////////

Loads the noise-added image and denoise it with a third-order FoE prior.
The experiment in this demo is described in the following paper:

Hiroshi Ishikawa, "Higher-Order Clique Reduction without Auxiliary Variables,"
In CVPR2014, Columbus, Ohio, June 23-28, 2014.

It also uses techniques described in the following papers:

Hiroshi Ishikawa, "Transformation of General Binary MRF Minimization to the First Order Case,"
IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 33, no. 6, pp. 1234-1249,
June 2011.

Hiroshi Ishikawa, "Higher-Order Clique Reduction in Binary Graph Cut,"
In CVPR2009, Miami Beach, Florida. June 20-25, 2009.

FoE is described in the following paper:

Stefan Roth and Michael J. Black, "Fields of Experts: A Framework for Learning Image Priors,"
In CVPR2005, San Diego, California, June 20-25, 2005, p.p. II:860-867.

This software requires the QPBO software by Vladimir Kolmogorov available
at http://pub.ist.ac.at/~vnk/software.html

This software has been tested on Windows 7 (x64) with Visual Studio 2010,
Ubuntu 12.04 with g++ 4.6.3, and Ubuntu 12.10 with g++ 4.8.1.
Any report on bugs and results of trying on other platforms is appreciated.

//////////////////////////////////////////////////////////////////////////////////////////////////*/

#include "ELC/ELC.h"
#include "QPBO/QPBO.h"
#include "Image.h"
#include <string>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <algorithm>
#include <iostream>

using namespace ELCReduce;

// Image filenames
std::string NOISE_ADDED = "Images2/noise"; // The given noisy image.
std::string ORIGINAL = "Images2/original"; // For PSNR calculation
std::string RESULT = "Images2/result";  // To save the denoised image.
std::string LOG = "log2/log";
FILE* Log;
/*
const char* NOISE_ADDED = "test001_020.pgm"; // The given noisy image.
const char* ORIGINAL = "test001.pgm"; // For PSNR calculation
const char* RESULT = "test001_result.pgm"; // To save the denoised image.
*/

// Parameters
typedef int REAL; // Type for the energy value (only integer has been tested)
const REAL sigma = 20;	// Sigma of the noise added
const int mult = 1000;	// Multiplier for representing real values by integers
const int Ecycle = 20;	// Number of iterations to check energy decrease.
const int maxIt = 300;	// Maximum number of iteration.
const double stopThreshold = 100.0; // Stop if the energy has changed less than this after Ecycle iterations.
const double switchThreshold = 900.0;
const int Bcycle = 30;	// Number of iterations to renew the blurred image.

// direction vectors
const int vx[4] = {0, 1, 0, 1};
const int vy[4] = {0, 0, 1, 1};

// FoE experts provided by Stefan Roth
double alpha[3] = {0.586612685392731, 1.157638405566669, 0.846059486257292};
double expert[3][4] = {
	{-0.0582774013402734, 0.0339010363051084, -0.0501593018104054, 0.0745568557931712},
	{0.0492112815304123, -0.0307820846538285, -0.123247230948424, 0.104812330861557},
	{0.0562633568728865, 0.0152832583489560, -0.0576215592718086, -0.0139673758425540}
};


template<typename T> T square(const T& t) {return t * t;}


double getFoELocal(unsigned int h[4])
{
	double e = 0;
	for (int j = 0; j < 3; j++)
	{
		double s = 0;
		for (int k = 0; k < 4; k++)
			s += expert[j][k] * h[k];
		e += alpha[j] * log(1 + 0.5 * square(s));
	}
	return e;
}


double unaryEnergy(int x, int data, int sigma)
{
	return (double)square(x - data) / (square(sigma) * 2);
}


template<typename F>
void addEnergy(F& f, const image& cur, const image& im, const image& proposal, int W, int H, int N, REAL mult, REAL sigma)
{
	for (int y = 0; y + 1 < H; y++)	// FoE prior
		for (int x = 0; x + 1 < W; x++)
		{
			int nds[4];
			for (int j = 0; j < 4; j++)
				nds[j] = x + vx[j] + (y + vy[j]) * W;
			REAL E[16];
			for (int i = 0; i < 16; i++) // i = 0000 (binary) means all current, i = 1000 means only node 0 is proposed,
			{                                      // i = 1001 means node 0 and 3 are proposed, etc.
				unsigned int h[4];
				int b = 8;
				for (int j = 0; j < 4; j++)
				{
					h[j] = (i & b) ? proposal.buf[nds[j]] : cur.buf[nds[j]];
					b >>= 1;
				}
				E[i] = (REAL)(getFoELocal(h) * mult);
			}
			f.AddHigherTerm(4, nds, E);
		}
	for (int j = 0; j < N; j++) // Data term
	{
		double e0 = unaryEnergy(cur.buf[j], im.buf[j], sigma);
		double e1 = unaryEnergy(proposal.buf[j], im.buf[j], sigma);
		f.AddUnaryTerm(j, (REAL)(e0 * mult), (REAL)(e1 * mult));
	}
}


void buildQPBF(QPBO<REAL>& qpbo, const image& cur, const image& im, const image& proposal, int W, int H, int N, REAL mult, REAL sigma, int mode)
{
    clock_t begin_add = clock();
	PBF<REAL> pbf(W * H * 10);
	addEnergy(pbf, cur, im, proposal, W, H, N, mult, sigma);
	printf("%.3f\t", (float)(clock() - begin_add) / CLOCKS_PER_SEC);
	fprintf(Log,"%.3f\t", (float)(clock() - begin_add) / CLOCKS_PER_SEC);

	if (mode == 0)
	{	// mode 0: reduce only ELCs, convert the rest with HOCR
        clock_t begin_reduceHigher = clock();
		pbf.reduceHigher(); // Reduce ELCs
        //pbf.printStar();
        printf("%.3f\t", (float)(clock() - begin_reduceHigher) / CLOCKS_PER_SEC);
        fprintf(Log,"%.3f\t", (float)(clock() - begin_reduceHigher) / CLOCKS_PER_SEC);

        clock_t begin_toQuadratic = clock();
		PBF<REAL> qpbf(W * H * 10);
		pbf.toQuadratic(qpbf, N); // Reduce what is left with HOCR. N is the ID for new variable.
		printf("%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);
        fprintf(Log,"%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);

		pbf.clear(); // free memory
		//clock_t begin_convert = clock();
		qpbf.convert(qpbo, N, Log); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
		//printf("%.4f\t", (float)(clock() - begin_convert) / CLOCKS_PER_SEC);
	}
	else if (mode == 1)
	{	// mode 1: reduce all higher-order terms using the approximation
		clock_t begin_reduceHigherApprox = clock();
		pbf.reduceHigherApprox();
		printf("%.3f\t\t", (float)(clock() - begin_reduceHigherApprox) / CLOCKS_PER_SEC);
		fprintf(Log,"%.3f\t\t", (float)(clock() - begin_reduceHigherApprox) / CLOCKS_PER_SEC);

		pbf.convert(qpbo, N, Log); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
	else if (mode == 2)
	{	// mode 2: use only HOCR
	    clock_t begin_toQuadratic = clock();
		PBF<REAL> qpbf(W * H * 10);
		pbf.toQuadratic(qpbf, N); // Reduce to Quadratic pseudo-Boolean function. N is the ID for new variable.
		printf("%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);
        fprintf(Log,"%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);

		pbf.clear(); // free memory
		qpbf.convert(qpbo, N, Log); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
	else if (mode == 3)
	{	// mode 3: use only Tin
	    clock_t begin_toQuadratic = clock();
		PBF<REAL> qpbf(W * H * 10);
		pbf.toQuadratic_Tin(qpbf, N, W); // Reduce to Quadratic pseudo-Boolean function. N is the ID for new variable.
		printf("%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);
        fprintf(Log,"%.3f\t", (float)(clock() - begin_toQuadratic) / CLOCKS_PER_SEC);

		pbf.clear(); // free memory
		qpbf.convert(qpbo, N, Log); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
}



double getEnergy(image& f, const image& data, int W, int H, REAL sigma)
{
	double E = 0;
	unsigned int h[4];
	for (int y = 0; y + 1 < H; y++)
		for (int x = 0; x + 1 < W; x++)
		{
			for (int j = 0; j < 4; j++)
				h[j] = f(x + vx[j], y + vy[j]);
			E += getFoELocal(h);
		}
	for (int j = 0; j < H * W; j++)
		E += unaryEnergy(f.buf[j], data.buf[j], sigma);
	return E;
}


double getPSNR(const image& d, const image& o, int N)
{
	int s = 0;
	for (int i = 0; i < N; i++)
		s += square(d.buf[i] - o.buf[i]);
	return 20 * log10(255.0 / sqrt((double)s / N));
}


double rand01()
{
	return ((double)rand()-0.00000000001) / RAND_MAX;
}

void create_images(int dim_size){
	int mu = 0, sigma = 20, step = 255/dim_size;

	int i, j, temp, noise_added;
	int width = dim_size, height = dim_size;

    // initialize random seed:
    //srand (time(NULL));
    // generate number between 1 and 10:
    // rand() % 10 + 1;

    std::default_random_engine generator;
    std::normal_distribution<double> gaussian(mu, sigma);
    //std::uniform_int_distribution<double> uniform(0, 99); // from 0 to 99 inclusive

	FILE* pgm_original;
	FILE* pgm_noise;
	pgm_original = fopen( (ORIGINAL + "_" + std::to_string(dim_size) + ".pgm").c_str(), "wb");
	pgm_noise    = fopen( (NOISE_ADDED + "_" + std::to_string(dim_size) + ".pgm").c_str(), "wb");

	// Writing Magic Number to the File
	fprintf(pgm_original, "P2\n\n");
	fprintf(pgm_noise, "P2\n\n");

	// Writing Width and Height
	fprintf(pgm_original, "%d %d\n", width, height);
    fprintf(pgm_noise, "%d %d\n", width, height);

	// Writing the maximum gray value
	fprintf(pgm_original, "255\n");
	fprintf(pgm_noise, "255\n");
	for (i = 0; i < height; i++) {
        /*
        if (i % 2 == 1)
            temp = 255 - 3*(i+1);
        else
            temp = 255 - 3*i;
        */
        //temp = 255 - step*i;
		for (j = 0; j < width; j++) {
			temp = 255 - std::min(i, std::min(j, std::min(height-i, width-j)));// * 2;
			// Writing the gray values in the 2D array to the file
			fprintf(pgm_original, "%d ", temp);

			// Add Gaussian Noise
			noise_added = std::max(0, std::min(255, temp + int(gaussian(generator)) ) );
			/*
			// Add Salt and Pepper Noise
			rand = uniform(generator);
			if (rand < 5)
				noise_added = 0;
			else if (rand > 94)
				noise_added = 255;
			else
				noise_added = temp;
			*/
			fprintf(pgm_noise, "%d ", noise_added);
		}
		fprintf(pgm_original, "\n");
		fprintf(pgm_noise, "\n");
	}
	fclose(pgm_original);
	fclose(pgm_noise);
}

void plot(int mode, int start, int end_)
{
	int k = 0;
	float avg[20][7];
	Log = fopen( (LOG + "_" + std::to_string(mode) + ".txt").c_str(), "w");
	for(int dim_size = start; dim_size <= end_; dim_size+=10)
    {
    	//if (dim_size == 210) dim_size = 300;
        std::ifstream s( (LOG + "_" + std::to_string(mode) + "_" + std::to_string(dim_size) + ".txt").c_str() );
		if(!s)
            printf("File %s is empty!\n", (LOG + "_" + std::to_string(mode) + "_" + std::to_string(dim_size) + ".txt").c_str());
        else{
            // Skip 4 lines
            char b[5000];
            for(int i=0;i<4;i++) s.getline(b, 4999);

            int counter = 0;
            float build, elc, hocr, solve, vars, quadratics, sum[] = {0,0,0,0,0,0};
            while(s){
                s >> build;
                if (mode!=2) s >>  elc;
				if (mode!=1) s >> hocr;
				for(int i=0;i<3;i++) s >> quadratics;
				s >> solve;
				for(int i=0;i<6;i++) s >> vars;
				s.getline(b, 4999);
                sum[0] += build;
                if (mode!=2) sum[1] += elc;
                if (mode!=1) sum[2] += hocr;
                sum[3] += solve;
                sum[4] += vars;
                sum[5] += quadratics;
                counter++;
                printf("%d, %d, %d, %d, %d\n", build, hocr, quadratics, solve, vars);
            }
            //last entry was added twice
            counter--;
            avg[k][0] = (sum[0] - build)/counter;
			if (mode!=2) avg[k][1] = (sum[1] -   elc)/counter;
			if (mode!=1) avg[k][2] = (sum[2] -  hocr)/counter;
            avg[k][3] = (sum[3] - solve)/counter;
            avg[k][4] = (sum[4] -  vars)/counter;
            avg[k][5] = dim_size * dim_size;
            avg[k][6] = (sum[5] -  quadratics)/counter;
            k++;
        }
    }
    ///////////////////////////////////////--BUILD--/////////////////////////////////////////

    fprintf( Log,"plot( ["); printf("plot( [");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.0f, ", avg[i][5]); printf("%.0f, ", avg[i][5]); }
    fprintf( Log,"%.0f], ", avg[k-1][5]); printf("%.0f], ", avg[k-1][5]);

    fprintf( Log,"["); printf("[");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][0]); printf("%.3f, ", avg[i][0]); }
    fprintf( Log,"%.3f] );\n", avg[k-1][0]); printf("%.3f] );\n", avg[k-1][0]);

    ///////////////////////////////////////---ELC---/////////////////////////////////////////
	if (mode!=2)
	{
		fprintf( Log,"plot( ["); printf("plot( [");
		for(int i=0;i<k-1;i++){ fprintf( Log,"%.0f, ", avg[i][5]); printf("%.0f, ", avg[i][5]); }
		fprintf( Log,"%.0f], ", avg[k-1][5]); printf("%.0f], ", avg[k-1][5]);

		fprintf( Log,"["); printf("[");
		for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][1]); printf("%.3f, ", avg[i][1]); }
		fprintf( Log,"%.3f] );\n", avg[k-1][1]); printf("%.3f] );\n", avg[k-1][1]);
	}
    ///////////////////////////////////////--HOCR--//////////////////////////////////////////
	if (mode!=1)
	{
		fprintf( Log,"plot( ["); printf("plot( [");
		for(int i=0;i<k-1;i++){ fprintf( Log,"%.0f, ", avg[i][5]); printf("%.0f, ", avg[i][5]); }
		fprintf( Log,"%.0f], ", avg[k-1][5]); printf("%.0f], ", avg[k-1][5]);

		fprintf( Log,"["); printf("[");
		for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][2]); printf("%.3f, ", avg[i][2]); }
		fprintf( Log,"%.3f] );\n", avg[k-1][2]); printf("%.3f] );\n", avg[k-1][2]);
	}
    ///////////////////////////////////////--SOLVE--/////////////////////////////////////////

    fprintf( Log,"\nplot( ["); printf("\nplot( [");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][4]); printf("%.3f, ", avg[i][4]); }
    fprintf( Log,"%.3f], ", avg[k-1][4]); printf("%.3f], ", avg[k-1][4]);

    fprintf( Log,"["); printf("[");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][3]); printf("%.3f, ", avg[i][3]); }
    fprintf( Log,"%.3f] );\n", avg[k-1][3]); printf("%.3f] );\n", avg[k-1][3]);

    ///////////////////////////////////////--Quadratics--/////////////////////////////////////////

    fprintf( Log,"\nplot( ["); printf("\nplot( [");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.0f, ", avg[i][4]); printf("%.0f, ", avg[i][4]); }
    fprintf( Log,"%.0f], ", avg[k-1][4]); printf("%.0f], ", avg[k-1][4]);

    fprintf( Log,"["); printf("[");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][6]); printf("%.3f, ", avg[i][6]); }
    fprintf( Log,"%.3f] );\n", avg[k-1][6]); printf("%.3f] );\n", avg[k-1][6]);

    ///////////////////////////////////////--#Variables--/////////////////////////////////////////

    fprintf( Log,"\nplot( ["); printf("\nplot( [");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.0f, ", avg[i][5]); printf("%.0f, ", avg[i][5]); }
    fprintf( Log,"%.0f], ", avg[k-1][5]); printf("%.0f], ", avg[k-1][5]);

    fprintf( Log,"["); printf("[");
    for(int i=0;i<k-1;i++){ fprintf( Log,"%.3f, ", avg[i][4]); printf("%.3f, ", avg[i][4]); }
    fprintf( Log,"%.3f] );\n", avg[k-1][4]); printf("%.3f] );\n", avg[k-1][4]);
}

const char* modename[] = {"ELC+HOCR", "Approx. ELC", "HOCR", "Tin"};
int main(int argc, char *argv[])
{
//	int random_seed;
//	std::cin >> random_seed;
//	srand(random_seed);

    int mode = 3;
	if (argc == 2)
		mode = atoi(argv[1]);
	double Erec[Ecycle];
	// mode 0: reduce only ELC, convert the rest with HOCR
	// mode 1: reduce all higher-order terms using the approximation
	// mode 2: use only HOCR

	if (false)
		plot(mode, 10, 100);
	else
	{
	for(int dim_size = 170; dim_size <= 170; dim_size+=10)
    {
        Log = fopen( (LOG + "_" + std::to_string(mode) + "_" + std::to_string(dim_size) + ".txt").c_str(), "w");

		printf (    "Mode\t%d (%s)\tOriginal image:\t%s\tNoise-added image:\t%s\n", mode, modename[mode], (ORIGINAL + "_" + std::to_string(dim_size) + ".pgm").c_str(), (NOISE_ADDED + "_" + std::to_string(dim_size) + ".pgm").c_str());
        fprintf(Log,"Mode\t%d (%s)\tOriginal image:\t%s\tNoise-added image:\t%s\n", mode, modename[mode], (ORIGINAL + "_" + std::to_string(dim_size) + ".pgm").c_str(), (NOISE_ADDED + "_" + std::to_string(dim_size) + ".pgm").c_str());

        //create_images(dim_size);
        image im, org, blr;
        im.readPGM( (NOISE_ADDED + "_" + std::to_string(dim_size) + ".pgm").c_str() );
        org.readPGM((ORIGINAL    + "_" + std::to_string(dim_size) + ".pgm").c_str() );

        if (im.empty() || org.empty())
        {
            printf("Error. Cannot load the image.\n");
            return 0;
        }
        int W = im.W, H = im.H, N = W * H;
        image cur(im);
        image proposal(W, H);

        printf("Initial Energy =\t%.1f\tInitial PSNR =\t%.3f\n", (double)getEnergy(im, im, W, H, sigma), getPSNR(im, org, N));
        fprintf(Log,"Initial Energy =\t%.1f\tInitial PSNR =\t%.3f\n", (double)getEnergy(im, im, W, H, sigma), getPSNR(im, org, N));
        if (mode == 0){
            printf("\nt_build\tt_ELC\tt_HOCR\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
            fprintf(Log,"\nt_build\tt_ELC\tt_HOCR\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
        }
        else if (mode == 1){
            printf("\nt_build\tt_ELCapprox\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
            fprintf(Log,"\nt_build\tt_ELCapprox\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
        }
        else if (mode == 2){
            printf("\nt_build\tt_HOCR\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
            fprintf(Log,"\nt_build\tt_HOCR\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
        }
        else if (mode == 3){
            printf("\nt_build\tt_Tin\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
            fprintf(Log,"\nt_build\tt_Tin\t\trange\t\t\t#quadratics\tt_solve\tIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
        }
        clock_t begin = clock(), begin_build, begin_solve;
        for (int i = 0; i < maxIt; i++)
        {
            if ((i % Bcycle) == 0)
                cur.gaussianblur(blr, .5625);
            if ((i % 2) == 0)
                for (int j = 0; j < N; j++)
                    proposal.buf[j] = std::min((REAL)255, std::max((REAL)0, (REAL)blr.buf[j] + (REAL)((rand01()-0.5) * sigma * 3)));
            else
                for (int j = 0; j < N; j++)
                    proposal.buf[j] = (int)(rand01() * 256);

            QPBO<REAL> qpbo(N * 10, N * 20);
            //begin_build = clock();
            buildQPBF(qpbo, cur, im, proposal, W, H, N, mult, sigma, mode);
            //printf("%.3f\t", (float)(clock() - begin_build) / CLOCKS_PER_SEC);

            begin_solve = clock();
            qpbo.MergeParallelEdges();
            qpbo.Solve();
            qpbo.ComputeWeakPersistencies();
            printf("%.3f\t", (float)(clock() - begin_solve) / CLOCKS_PER_SEC);
            fprintf(Log,"%.3f\t", (float)(clock() - begin_solve) / CLOCKS_PER_SEC);

            int labeled = 0;
            for (int j = 0; j < N; j++)
            {
                int res = qpbo.GetLabel(j);
                if (res == 1)
                    cur.buf[j] = proposal.buf[j];
                if (res >= 0)
                    labeled++;
            }
            double psnr = getPSNR(cur, org, N);
            double te = getEnergy(cur, im, W, H, sigma);
            int elapsed = clock() - begin;
            printf("%d\t%.3f\t", i, (float)elapsed / CLOCKS_PER_SEC);
            printf("%.2f\t", labeled * 100.0 / N);
            printf("%.1f\t%.3f\t%d\n", te, psnr, qpbo.GetNodeNum());
            fprintf(Log,"%d\t%.3f\t", i, (float)elapsed / CLOCKS_PER_SEC);
            fprintf(Log,"%.2f\t", labeled * 100.0 / N);
            fprintf(Log,"%.1f\t%.3f\t%d\n", te, psnr, qpbo.GetNodeNum());
            fflush(stdout);
            //if (i >= Ecycle && Erec[i % Ecycle] - te < switchThreshold)
            //	mode = 0;
            if (i >= Ecycle && Erec[i % Ecycle] - te < stopThreshold)
                break;
            Erec[i % Ecycle] = te;
        }
        const char* out = (RESULT + "_" + std::to_string(mode) + "_" + std::to_string(dim_size) + ".pgm").c_str();
        cur.writePGM(out);
        printf("Denoised image written to\t%s\n", out);
        fprintf(Log,"Denoised image written to\t%s\n", out);
        fclose(Log);
    }
	}
    return 0;
}

/*
int main(int argc, char *argv[])
{
    int W = 3, H = 3, N = W * H;
    image im(W, H), proposal(W, H);

    im.buf = {201, 180, 195, 183, 209, 210, 173, 171, 162};
    proposal.buf = {242, 230, 255, 234, 252, 239, 238, 204, 241};

    image cur(im);

    PBF<REAL> pbf(W * H * 10);
	addEnergy(pbf, cur, im, proposal, W, H, N, mult, sigma);
	pbf.shrink();
	pbf.print();
}
*/

