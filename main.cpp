#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>


/*
SORUNLAR
1)
*/

using namespace std;

//float min = 1000;
float minimum = 1000;
float maximum = 0;



const int pop_Size = 10;
const float KE_Loss_Rate = 0.2;
const float MoleColl = 0.4;
const int buffer = 0;
int Initial_KE = 1000;
const int alpha = 5;
const int beta = 15000;

const int num_node = 402;
const int num_molecule = 10;
const int num_parents = 5;
const int range = 800;
const int num_facility = 6;
const int iter = 5;





const float mutation_rate = 0.3;
const float replacement_prob = 0.3;
const int num_iter = 20;
const float prob_repairing = 0.5;
string problem_demand = "sjc402_demand.txt";
string problem_distance = "sjc402_distance.txt";
string problem_result = "sjc402-1p-max_imp.txt";
string improving="a";
string xover_operator="b";

struct node {
	double x_cord;
	double y_cord;
	int demand;
	int satisfied_demand;
};

node Array_node[num_node];
float distance_node[num_node][num_node];

struct molecule {
	int soln[num_node];
	int satisfied[num_node];
	int PE;///Potential Energy = fitness value
	int KE;///Kinetic Energy
	int NumHit;
	int MinStruct[num_node];
	int MinPE;
	int MinHit;
	float prob;
	int index;
};

struct population {
	molecule chros[num_molecule];
};

struct parent {
	///molecule parents[num_parents];
	vector<molecule> parents[num_parents];
};

parent mating_pool;

population initial;

void init_pop(void)  ///Generating Initial Population
{
	for (int i = 0; i < num_molecule; i++)
    {
		int counter = 0;
		while (counter <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (initial.chros[i].soln[rand_num] == 0)
            {
				initial.chros[i].soln[rand_num] = 1;
				counter += 1;
			}
		}
		initial.chros[i].index = i;
        initial.chros[i].KE = Initial_KE;
        initial.chros[i].NumHit = 0;
        initial.chros[i].MinHit = 0;
	}
}

int total_demand[num_node];

void demand_satisfied_node() { ///Calculates each nodes' satisfaction portion in case facility opened

	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < num_node; j++)
		{
			if (distance_node[i][j] < range) {
				total_demand[i] += Array_node[j].demand;
			}
		}
		Array_node[i].satisfied_demand = total_demand[i] + Array_node[i].demand;
	}
}

void ordering() { ///Ordering nodes according to their demand_satisfied_node()
	int i, j;
	node temp;

	for (i = 0; i < num_node; i++)
	{
		for (j = i + 1; j < num_node; j++)
		{
			if (Array_node[i].satisfied_demand < Array_node[j].satisfied_demand)
			{
				temp = Array_node[i];
				Array_node[i] = Array_node[j];
				Array_node[j] = temp;
				i = 0;
			}
		}
	}
}

void satisfied_vec_func(molecule &x)///satisfied vector calculation
{
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}

	for (int i = 0; i < num_node; i++)
    {
		for (int j = 0; j < num_node; j++)
        {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range))
            {
				x.satisfied[j] = 1;
			}
			else
            {
				if(x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

void satisfied_vec_func_offspring(molecule &x) { //satisfied vector calculation for offsprings
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

int fitness_func(molecule x) { //Fitness function calculation
	x.PE = 0;
	int sum = 0;
	for (int i = 0; i < num_node; i++) {
		sum = x.satisfied[i] * Array_node[i].demand + sum;
		/*
		if(x.satisfied[i]==1)
        {
            cout << Array_node[i].demand << " ";///kon kon demand value add kore fitness value passi seta
        }
        */
	}
	///cout << "Sum=";
	///cout << sum << endl;
	x.PE = sum;
	return x.PE;
}

void prob_func(population &gen) /// Calculating probabilities of chromosomes ***Note: Min Fitness Prob=0***
{
	for (int i = 0; i < num_molecule; i++)
    {
		if (gen.chros[i].PE < minimum)
        {
			minimum = gen.chros[i].PE;
		}
		else if (gen.chros[i].PE > maximum)
		{
			maximum = gen.chros[i].PE;
		}
	}
	///cout << maximum << endl;
	///cout << minimum << endl;
	for (int i = 0; i < num_molecule; i++) {
		gen.chros[i].prob = (gen.chros[i].PE - minimum) / (maximum - minimum);
	}
	float sum = 0;
	for (int i = 0; i < num_molecule; i++)
	{
		sum += gen.chros[i].prob;
	}
	///cout << sum;
	gen.chros[0].prob = gen.chros[0].prob / sum;
	for (int i = 1; i < num_molecule; i++)
	{
		gen.chros[i].prob = gen.chros[i - 1].prob + (gen.chros[i].prob / sum); //Finding Cumulative Probabilities
	}
}

float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

/*
void parent_selection(population gen) { ///Selecting # of distinct parents out of population

	int m = -1;
	int iteration = 0;

	for (int h = 0; h < num_parents + iteration; h++)
	{
		m++;
		float rnd = RandomFloat(0, 1);
		float diff = 1;
		for (int j = 0; j < num_molecule; j++)
		{
			if (((gen.chros[j].prob - rnd) > 0) && ((gen.chros[j].prob - rnd) < diff)) {
				diff = abs(gen.chros[j].prob - rnd);
                mating_pool.parents[m] = gen.chros[j];
                ///mating_pool.parents.push_back(gen.chros[j]);
				mating_pool.parents[m].index = gen.chros[j].index;
				///mating_pool.parents[m].index = gen.chros[j].index;


				for (int k = 0; k < m; k++)
				{
					if (mating_pool.parents[k].index == mating_pool.parents[m].index)
					{
						m = m - 1;
						iteration++;
					}
				}
			}
		}
	}
}

*/

molecule crossover(molecule parent1, molecule parent2) //Crossover Operator(1-point)
{
	xover_operator = "1-Point Crossover";
	int k = 188; //Cut Point
	molecule offspring;
	for (int i = 0; i < k; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule crossover_2p(molecule parent1, molecule parent2) //Crossover Operator(2-point)
{
	xover_operator = "2-Point Crossover";
	int k1 = 180; //Cut Point
	int k2 = 246;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i = k1; i < k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule tri_crossover(molecule parent1, molecule parent2, molecule parent3) //3-Parent Crossover
{
	int k1 = 182;
	int k2 = 265;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i =k1; i <k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent3.soln[i];
	}
	return offspring;
}

molecule satisfied_new(molecule x) { //satisfied vector calculation
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
	return x;
}

molecule max_imp_repairing(molecule x) {
	improving = "Maximum Improving";
	int sum = 0;
	int max_demand = 0;
	int min_demand = 10000;
	int index = 0;
	x = satisfied_new(x);
	x.PE = fitness_func(x);
	int current_fitness = x.PE;
	molecule y;
	y = x;
	int deviation[num_node];

	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}
	if (sum < num_facility) {
		while (sum < num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 0)
				{
					y.soln[i] = 1;
					y=satisfied_new(y);
					y.PE=fitness_func(y);
					deviation[i] = y.PE - current_fitness;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if (deviation[i] > max_demand)
				{
					max_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 1;
			sum++;
		}
	}

	if (sum > num_facility) {
		while (sum > num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 1)
				{
					y.soln[i] = 0;
					y=satisfied_new(y);
					y.PE = fitness_func(y);
					deviation[i] = current_fitness-y.PE;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if ((deviation[i] <min_demand)&&(deviation[i]>0))
				{
					min_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 0;
			sum--;
		}
	}
	return x;
}

molecule repairing(molecule x) { //Randomly Repairing
	improving = "Random Repairing";
	int sum = 0;
	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}

	if (sum < num_facility) {
		while (sum <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 0) {
				x.soln[rand_num] = 1;
				sum += 1;
			}
		}
	}
	if (sum > num_facility) {
		while (sum >= num_facility + 1) {
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 1) {
				x.soln[rand_num] = 0;
				sum -= 1;
			}
		}
	}
	return x;
}

molecule alt_repairing(molecule x) { //Repairing according to demand_satisfied_node
	improving = "Greedy Repairing";
	int sum = 0;
	float rand_num = 0;
	for (int i = 0; i < num_node; i++)
	{

		sum += x.soln[i];
	}

	int maks = 0;

	while (sum < num_facility) {
		for (int i = 0; i < num_node; i++)
		{
			rand_num = RandomFloat(1, 0);
			if ((x.soln[i] == 0) && (rand_num < prob_repairing))
			{
				x.soln[i] = 1;
				sum++;
			}
		}
	}

	while (sum > num_facility)
	{
		for (int i = num_node; i > -1; i--)
		{
			rand_num = RandomFloat(1, 0);

			if ((x.soln[i] == 1) && (rand_num < prob_repairing))
			{
				x.soln[i] = 0;
				sum--;
			}
		}
	}

	return x;
}

molecule mutation(molecule x) { ///Mutation operator (swap two elements of array)
	float rand_num = RandomFloat(1, 0);
	int rnd = rand() % num_node;
	if (rand_num < mutation_rate) {
		for (int i = 0; i < num_node; i++)
		{
			swap(x.soln[rnd], x.soln[rand() % num_node]);
		}
	}
	return x;
}

int main()
{
	int bekle;
	int mak_out = 0;
	ofstream result("Result.txt");
	double timer=0;
	double timer_out = 0;
	int sum_time = 0;
	int sum_best = 0;


		int start_s = clock();
		// Start timer

		srand((int)time(NULL));

		ofstream out("sjc402-1p-max_imp.txt");
		ifstream get_demand; //Get demand values from problem set

		get_demand.open("sjc402_demand.txt");
		for (int i = 0; i < num_node; i++)
		{
			get_demand >> Array_node[i].demand;
		}

		get_demand.close();

		ifstream get_distance; //Get distance values from problem set

		get_distance.open("sjc402_distance.txt");

		for (int i = 0; i < num_node; i++)
		{
			for (int j = 0; j < num_node; j++)
			{
				get_distance >> distance_node[i][j];
			}
		}

		get_distance.close();

		demand_satisfied_node(); //Calculate each node's portion of satisfaction

		//ordering(); //Order nodes according to portion of satisfaction

		init_pop(); //Initialize a population



        for (int i = 0; i < num_molecule; i++) //Calculating satisfied vector
        {
            satisfied_vec_func(initial.chros[i]);
        }

        ///Initialize MinPE
        for (int j = 0; j < num_molecule; j++) //Calculate fitness values of population
        {
            initial.chros[j].PE = fitness_func(initial.chros[j]);
            initial.chros[j].MinPE = initial.chros[j].PE;
        }

        ///Initialize MinStruct
        for (int i = 0; i < num_molecule; i++)
        {
            for (int j = 0; j < num_node; j++)
            {
                initial.chros[i].MinStruct[j] = initial.chros[i].soln[j];
            }
        }

        prob_func(initial); /// Calculate Probabilities for mating


        int iter = 0;
        cout << endl;
        cout << iter << "th Generation" << endl;
        cout << endl;
        cout << "\t" "Chromosome" << "\t" << "Fitness" << "\t" << "Probabilities" << endl << endl;

        for (int i = 0; i < num_molecule; i++)
        {
            cout << initial.chros[i].index << "\t";
            for (int j = 0; j < num_node; j++)
            {
                ///cout << initial.chros[i].MinStruct[j];
                cout << initial.chros[i].satisfied[j];
            }
            cout << "\t";
            cout << initial.chros[i].PE << "\t";
            cout << initial.chros[i].prob << "\t";
            cout << initial.chros[i].KE << "\t";
            cout << initial.chros[i].NumHit << "\t";
            cout << initial.chros[i].MinHit << "\t";
            cout << initial.chros[i].MinPE << "\t";
            cout << endl;
        }

        ///float b = RandomFloat(0,1);
        ///cout << b;

        /*
        parent_selection(initial);

			cout << "Selected Parent Indices: ";

			for (int i = 0; i < num_parents; i++)
			{
				cout << mating_pool.parents[i].index << "-";
			}
			cout << endl;

         */

        ///Iterations
        int i=0;

        while (i < iter)
        {
            float b = RandomFloat(0,1);
            if (b > MoleColl)
            {

            }
        }

        goto ende;


	ende: return 0;
}
