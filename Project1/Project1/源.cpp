#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
#define M 1500
#define K 40
#define V 11980
#define MaxRand 32761
#define MaxLine 200000
using namespace std;

int nmk[M][K], nkt[K][V], nktS[K];
double alpha = 1.25, beta = 0.01;
double phi[K][V], theta[M][K];
vector<int> Docment[M];
vector<int> Z[M];
map<string, int> Str2Int;


int iteration = 500;
double llhw(){
	double lgamma_alpha = lgamma(alpha);
	double lgamma_beta = lgamma(beta);



	double res = 0;
	res += K* lgamma(beta*V);
	for (int k = 0; k < K;++k)
		res -= lgamma(beta*V + nktS[k]);

	

	for (int t = 0; t < V; ++t){
		for (int k = 0; k < K; ++k){
			if (nkt[k][t] > 0)
				res += lgamma(beta + nkt[k][t]) - lgamma_beta;
		}
	}
	for (int m = 0; m < M; ++m)
	{
		res += (lgamma(alpha*K) - lgamma(alpha*K + Docment[m].size()));
		for (int k = 0; k < K; ++k){
			if (nmk[m][k] > 0)
				res += lgamma(alpha + nmk[m][k]) - lgamma_alpha;
		
		}

	}
	return res;
}
int  GetTotalNum(){
	int res = 0;
	for (int i = 0; i < M; ++i)
	{
		res += Docment[i].size();
	}
	return res;

}
void GetDocment(string path){
	ifstream in(path, ios::in);//以输入方式打开文件 
	int item;
	char line[MaxLine];
	int docidx = 0;
	int ItemInt = 0;
	while (!in.eof())
	{
		memset(line, 0, sizeof(line));

		in.getline(line,sizeof(line));
		int len = strlen(line);
		string item = "";
		istringstream buffin(line);
		
		vector<string> vec;

		while (!buffin.eof())
		{
			buffin >> item;
			vec.push_back(item);
		}
		Docment[docidx].resize(vec.size());
		int idx = 0;
		for (auto it : vec ){
			if (!Str2Int[it])
				Str2Int[it] = ItemInt++;
			Docment[docidx][idx++] = Str2Int[it];
		}
		docidx++;
		cout << docidx << endl;
	}
	cout << ItemInt << endl;
}

void Init(){

	memset(nmk, 0, sizeof(nmk));
	memset(nkt, 0, sizeof(nkt));
	memset(nktS, 0, sizeof(nktS));

	for (int i = 0; i < M; ++i)
	{
		Z[i].resize(Docment[i].size());

		for (int j = 0; j < Docment[i].size(); ++j)
		{
			int topic = rand() % K;

			Z[i][j] = topic;
			nmk[i][topic] ++;
			nkt[topic][Docment[i][j]]++;
			nktS[topic]++;
		}
	}
	
	cout << "init end" << endl;
		return;

}

void SetPhiAndTheta(){
	
	int nmkS[M];
	memset(nmkS, 0, sizeof(nmkS));
	for (int m = 0; m < M; ++m)
	{
		for (int k = 0; k < K; ++k){
			nmkS[m] += nmk[m][k];
		
		}
	
	}


	for (int m = 0; m < M; ++m){
		for (int k = 0; k < K; ++k){
			theta[m][k] = (nmk[m][k] + alpha) / (nmkS[m] + alpha*K);
			
		}
	}
	
	for (int k = 0; k < K; ++k){
	
		for (int t = 0; t < V; ++t){
			phi[k][t] = (nkt[k][t] + beta) / (nktS[k] + beta * V);
		}
	}



}

void OneIteration(){

	for (int i = 0; i < M; ++i){
		for (int j = 0; j < Docment[i].size(); ++j){
			int oldtopic = Z[i][j];
			nmk[i][oldtopic]--;
			nkt[oldtopic][Docment[i][j]]--;
			nktS[oldtopic]--;


			double prob[K];
			for (int k = 0; k < K; ++k){
				prob[k] = (nkt[k][Docment[i][j]] + beta) *(nmk[i][k] + alpha) / (nktS[k] + V * beta) ;
			}
	
			for (int k = 1; k < K; ++k)
				prob[k] += prob[k - 1];

			int newtopic;
			double ran = rand() % MaxRand;
			double normal =  ran / MaxRand * prob[K-1];
			

			for (int k = 0; k < K; ++k){
				if (prob[k] > normal )
				{
					newtopic = k;
					break;
				}
	
			}
			
			Z[i][j] = newtopic;
			nmk[i][newtopic]++;
			nkt[newtopic][Docment[i][j]]++;
			nktS[newtopic]++;

		}
	
	
	}
}

double  Compute_Perplexity(){
	double  result = 0;
	for (int m = 0; m < M; ++m)
	{
		for (int j = 0; j < Docment[m].size(); ++j)
		{
			int wi = Docment[m][j];
			double p[K];
			for (int i = 0; i < K; ++i)
				p[i] = 0;
			for (int k = 0; k < K; ++k)
				p[k] += theta[m][k] * phi[k][wi];
			double tepsum = 0;
			for (int k = 0; k < K; ++k)
				tepsum += p[k];

			result += log(tepsum);
		}
	}
	int TotalNum = GetTotalNum();
	result = (-result / TotalNum);

	return exp(result);
}
void Train(){
	SetPhiAndTheta();
	double perplexity = Compute_Perplexity();
	cout << "perplexity  " << perplexity << endl;
	for (int iter = 0; iter < iteration; ++iter){
		cout << "iter  " << iter << endl;
		OneIteration();
		double result = llhw();
		cout << result << endl;  
		SetPhiAndTheta();
		double perplexity = Compute_Perplexity();
		cout << "perplexity  " << perplexity << endl;
	}

	return;

}
void PrintTheta(string path){
	fstream out(path, ios::out);
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < K; ++j)
			out << theta[i][j] << " ";
		out << endl;
	}
}

void PrintResult(string path){
	ofstream out(path);
	for (int i = 0; i < M; ++i)
	{
		int Maxindex = 0; double MaxValue = theta[i][0];

		for (int k = 0; k < K; ++k)
		{
			if (theta[i][k] > MaxValue){
				MaxValue = theta[i][k];
				Maxindex = k;
			}

		}
		out << Maxindex << endl;

	}
	out.close();
}
int main(){
	int t = 20;
	srand(t);
	GetDocment("nips.train.txt");
	Init();
	Train();
	SetPhiAndTheta();

	PrintTheta("theta.txt");
	PrintResult("result1.txt");
	return 0;

}