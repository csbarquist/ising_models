#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define GRID_SIZE 10
#define TEMPERATURE 100

typedef struct Grid2D{
	int size;
	int position[GRID_SIZE][GRID_SIZE];
} Grid2D;

void init_grid(Grid2D* pg){
	int rand_int;  // storing random state
	
	// Set size for accessing later
	pg->size = GRID_SIZE;
	
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

int main(){
	//Init rng
	srand(time(NULL)); // seeds rand func
	
	// Init grid
	Grid2D grid;
	init_grid(&grid);
	print_grid(&grid);
	return 0;
}
