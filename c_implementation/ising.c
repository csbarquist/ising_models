#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EXTERNAL_FIELD 0       // strength of external magnetic field
#define INTERACTION_STRENGTH 1 // + for ferromagnetic, - for anti-feromagnetic
#define GRID_SIZE 40
#define BOLTZMANN_SCALE 10000

#define START_TEMP 1000
#define END_TEMP 0
#define TEMP_STEP -10
#define LOOP_COUNT 10000 // Number of updates per temp

typedef struct Grid2D {
  size_t size;
  int position[GRID_SIZE][GRID_SIZE];
  double field_strength;
  double interaction_strength;
} Grid2D;

typedef struct MeanStats {
  float m;
  float m2;
  float E;
  float E2;
} MeanStats;

typedef struct Deltas {
  int delta_E;
  int delta_M;
} Deltas;
void init_grid(Grid2D *pg, size_t grid_size, double field_strength,
               double interaction_strength) {
  int rand_int; //   storing random state

  // Set size for accessing later
  pg->size = grid_size;
  pg->field_strength = field_strength;
  pg->interaction_strength = interaction_strength;

  // Loop over rows and columns assiging 1 or 0 randomly
  // This is setting state of magnet high T random
  for (int row = 0; row < pg->size; row++) {
    for (int col = 0; col < pg->size; col++) {
      // randomly generate 1 or 0 and "double" it using bit shift
      // 1 - {2, 0} = {1, -1}
      rand_int = 1 - ((rand() % 2) << 1);
      pg->position[row][col] = rand_int;
    }
  }
  return;
}

void clear_screan() { printf("\e[1;1H\e[2J"); }

void print_grid(Grid2D *pg) {
  for (int row = 0; row < pg->size; row++) {
    for (int col = 0; col < pg->size; col++) {
      if (pg->position[row][col] < 0) {
        printf("x ");
      } else {
        printf("o ");
      }
    }
    printf("\n");
  }
}

int calc_magnetization(Grid2D *pg) {
  int mag = 0; //   init mag as 0 and sum
  for (int row = 0; row < pg->size; row++) {
    for (int col = 0; col < pg->size; col++) {
      mag += pg->position[row][col];
    }
  }
  return mag;
}

double site_energy_periodic(Grid2D *pg, int posy, int posx) {
  // Assume periodic boundary conditions, i.e., wrap
  // around when calculating energy at the edges.
  int site_magnetization = pg->position[posy][posx];

  double field_energy = -1 * pg->field_strength * site_magnetization;
  double interaction_energy = 0;

  int max_index = pg->size - 1;

  // x-axis interaction energy and edge case handling
  if (posx == 0) {
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy][max_index] * site_magnetization;
    interaction_energy -=
        pg->interaction_strength * pg->position[posy][1] * site_magnetization;
  } else if (posx == max_index) {
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy][max_index - 1] *
                          site_magnetization;
    interaction_energy -=
        pg->interaction_strength * pg->position[posy][0] * site_magnetization;
  } else {
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy][posx + 1] * site_magnetization;
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy][posx - 1] * site_magnetization;
  }
  // y-axis interaction energey and edge case handling
  if (posy == 0) {
    interaction_energy -= pg->interaction_strength *
                          pg->position[max_index][posx] * site_magnetization;
    interaction_energy -=
        pg->interaction_strength * pg->position[1][posx] * site_magnetization;
  } else if (posy == max_index) {
    interaction_energy -= pg->interaction_strength *
                          pg->position[max_index - 1][posx] *
                          site_magnetization;
    interaction_energy -=
        pg->interaction_strength * pg->position[0][posx] * site_magnetization;
  } else {
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy + 1][posx] * site_magnetization;
    interaction_energy -= pg->interaction_strength *
                          pg->position[posy][posx - 1] * site_magnetization;
  }

  return interaction_energy + field_energy;
}

double lattice_energy(Grid2D *pg) {
  // Calculates the total energy of the lattice using
  // E = -J*Sum(Si*Sj) - h*Sum(Si)

  double E = 0;
  for (int row = 0; row < pg->size; row++) {
    for (int col = 0; col < pg->size; col++) {
      E += site_energy_periodic(pg, row, col);
    }
  }
  return E;
}

bool boltzmann_coin(int delta_E, int temp) {
  // if return true then keep state
  int roll = rand() % BOLTZMANN_SCALE;
  int pass_val = floor(BOLTZMANN_SCALE * exp((-1.0 * delta_E) / temp));

  return roll < pass_val;
}

void flip_spin(Grid2D *pg, int row, int col) { pg->position[row][col] *= -1; }

Deltas *update_lattice(Grid2D *pg, Deltas *d, int temp) {
  int rand_row = rand() % GRID_SIZE;
  int rand_col = rand() % GRID_SIZE;
  int cur_E = site_energy_periodic(pg, rand_row, rand_col);
  int delta_E;

  flip_spin(pg, rand_row, rand_col); // Filp randomly chosen spin
  delta_E = site_energy_periodic(pg, rand_row, rand_col) - cur_E;

  // Update Diffs
  d->delta_E = delta_E;
  d->delta_M = 2 * pg->position[rand_row][rand_col];

  if (delta_E > 0 && !boltzmann_coin(delta_E, temp)) {
    flip_spin(pg, rand_row, rand_col);
    d->delta_E = 0;
    d->delta_M = 0;
  }

  return d;
}

void temp_loop_simulation(Grid2D *pg) {
  for (int temp = START_TEMP; temp) {
  }
}

int main() {
  // Init rng
  srand(time(NULL)); //   seeds rand func

  // Init grid
  Grid2D grid;
  init_grid(&grid, GRID_SIZE, EXTERNAL_FIELD, INTERACTION_STRENGTH);
  print_grid(&grid);

  double E = lattice_energy(&grid);
  printf("Lattice Energy %f\n", E);

  return 0;
}
