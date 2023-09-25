#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define EXTERNAL_FIELD 0  // strength of external magnetic field
#define INTERACTION_STRENGTH 1  // + for ferromagnetic, - for anti-feromagnetic
#define GRID_SIZE 10
#define TEMPERATURE 100


typedef struct Grid2D{
	size_t size;
	int position[GRID_SIZE][GRID_SIZE];
	double field_strength;
	double interaction_strength;
} Grid2D;


void init_grid(
		Grid2D* pg,
		size_t grid_size,
		double field_strength,
		double interaction_strength
	){
	int rand_int;  // storing random state
	
	// Set size for accessing later
	pg->size = grid_size;
	pg->field_strength = field_strength;
	pg->interaction_strength = interaction_strength;

	
	// Loop over rows and columns assiging 1 or 0 randomly
	// This is setting state of magnet high T random
	for(int row = 0; row < pg->size; row++){
		for(int col = 0; col < pg->size; col++){
			
			// randomly generate 1 or 0 and "double" it using bit shift
			// 1 - {2, 0} = {1, -1}
			rand_int = 1 - ((rand() % 2) << 1); 
			pg->position[row][col] = rand_int;	
		}
	}
	return;
}

void print_grid(Grid2D* pg){
	for(int row = 0; row < pg->size; row++){
		for(int col = 0; col < pg->size; col++){
        		printf("%d ", pg->position[row][col]);
		}
		printf("\n");
	} 
}

int calc_magnetization(Grid2D* pg){
	int mag = 0;  // init mag as 0 and sum
	for(int row = 0; row < pg->size; row++){
		for(int col = 0; col < pg->size; col++){
			mag += pg->position[row][col]; 
		}
	}
	return mag;
}

double site_energy_periodic(Grid2D* pg, int posy, int posx){
	// Assume periodic boundary conditions, i.e., wrap
	// around when calculating energy at the edges.
	int site_magnetization = pg->position[posy][posx];

	double field_energy = -1*pg->field_strength*site_magnetization;
	double interaction_energy = 0;

	int max_index = pg->size-1;

	// x-axis interaction energy and edge case handling
	if(posx==0){
		
		interaction_energy -= pg->interaction_strength*pg->position[posy][max_index]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[posy][1]*site_magnetization;
	}
	else if(posx==max_index){
		interaction_energy -= pg->interaction_strength*pg->position[posy][max_index-1]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[posy][0]*site_magnetization;
	}
	else{
		interaction_energy -= pg->interaction_strength*pg->position[posy][posx+1]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[posy][posx-1]*site_magnetization;
	}
	// y-axis interaction energey and edge case handling
	if(posy==0){
		
		interaction_energy -= pg->interaction_strength*pg->position[max_index][posx]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[1][posx]*site_magnetization;
	}
	else if(posy==max_index){
		interaction_energy -= pg->interaction_strength*pg->position[max_index-1][posx]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[0][posx]*site_magnetization;
	}
	else{
		interaction_energy -= pg->interaction_strength*pg->position[posy+1][posx]*site_magnetization;
		interaction_energy -= pg->interaction_strength*pg->position[posy][posx-1]*site_magnetization;
	}

	return interaction_energy + field_energy;
}

double lattice_energy(Grid2D* pg){
	// Calculates the total energy of the lattice using
	// E = -J*Sum(Si*Sj) - h*Sum(Si)

	double E = 0; 
	for(int row=0; row < pg->size; row++){
		for(int col=0; col < pg->size; col++){
			E += site_energy_periodic(pg, row, col);
		}
	}

	return E;
}

int main(){
	//Init rng
	srand(time(NULL)); // seeds rand func
	
	// Init grid
	Grid2D grid;
	init_grid(
		&grid,
		GRID_SIZE,
		EXTERNAL_FIELD,
		INTERACTION_STRENGTH
	);
	print_grid(&grid);

	double E = lattice_energy(&grid);
	printf("Lattice Energy %f\n", E);

	return 0;
}
