#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define OUT_FILE "./ising_out.csv"

#define EXTERNAL_FIELD 0       // strength of external magnetic field
#define INTERACTION_STRENGTH 1 // + for ferromagnetic, - for anti-feromagnetic
#define GRID_SIZE 40
#define BOLTZMANN_SCALE 10000

#define START_TEMP 3
#define END_TEMP 0
#define TEMP_STEP 0.05
#define LOOP_COUNT 1000 // Number of updates per temp

#define PRINT_FLAG true
#define PRINT_AFTER 500
#define PRINT_TIMER 1 // Seconds waited after printing

typedef struct Grid2D {
  size_t size;
  int position[GRID_SIZE][GRID_SIZE];
  double field_strength;
  double interaction_strength;
} Grid2D;

typedef struct MeanStats {
  double m;
  double m2;
  double E;
  double E2;
  int n; // Number of updates
} MeanStats;

typedef struct CurStats {
  double m;
  double E;
} CurStats;

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

void print_stats(MeanStats *ms, CurStats *cs, double temp) {
  printf("TEMP: %f\n<E>: %f\tE: %f\n<m>: %f\tm: %f\n", temp, ms->E, cs->E,
         ms->m, cs->m);
}

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

void init_mean_stats(MeanStats *ms, Grid2D *pg) {
  ms->m = (double)calc_magnetization(pg);
  ms->m2 = pow(ms->m, 2);
  ms->E = (double)lattice_energy(pg);
  ms->E2 = pow(ms->E, 2);
  ms->n = 1;
}

void init_cur_stats(CurStats *cs, Grid2D *pg) {
  cs->m = (double)calc_magnetization(pg);
  cs->E = (double)lattice_energy(pg);
}

bool boltzmann_coin(int delta_E, double temp) {
  // if return true then keep state
  int roll = rand() % BOLTZMANN_SCALE;
  int pass_val = floor(BOLTZMANN_SCALE * exp((-1.0 * delta_E) / temp));

  return roll < pass_val;
}

void flip_spin(Grid2D *pg, int row, int col) { pg->position[row][col] *= -1; }

void update_lattice(Grid2D *pg, CurStats *cs, double temp) {
  int rand_row = rand() % GRID_SIZE;
  int rand_col = rand() % GRID_SIZE;
  double cur_E = site_energy_periodic(pg, rand_row, rand_col);
  double delta_E;
  double delta_m;

  flip_spin(pg, rand_row, rand_col); // Filp randomly chosen spin
  delta_E = site_energy_periodic(pg, rand_row, rand_col) - cur_E;

  // Update Diffs
  delta_E = delta_E;
  delta_m = 2.0 * pg->position[rand_row][rand_col];

  // Heart of Metropolis method
  // If change in energey is < 0 automatically accept change
  // otherwise roll a weighted die.
  if (delta_E > 0 && !boltzmann_coin(delta_E, temp)) {
    flip_spin(pg, rand_row, rand_col);
    delta_E = 0;
    delta_m = 0;
  }

  cs->E += delta_E;
  cs->m += delta_m;
}

void update_stats(MeanStats *ms, CurStats *cs) {
  ms->m = ms->m * (ms->n) / (ms->n + 1) + cs->m / (ms->n + 1);
  ms->E = ms->E * (ms->n) / (ms->n + 1) + cs->E / (ms->n + 1);
  ms->m2 = ms->m2 * (ms->n) / (ms->n + 1) + pow(cs->m, 2) / (ms->n + 1);
  ms->E2 = ms->E2 * (ms->n) / (ms->n + 1) + pow(cs->E, 2) / (ms->n + 1);
}

void simulate_temp(Grid2D *pg, MeanStats *ms, CurStats *cs, double temp) {
  for (int i = 0; i < LOOP_COUNT; i++) {
    update_lattice(pg, cs, temp);
    update_stats(ms, cs);

    if (PRINT_FLAG && (i % PRINT_AFTER == 0)) {
      clear_screan();
      print_grid(pg);
      print_stats(ms, cs, temp);
      sleep(PRINT_TIMER);
    }
  }
}

void init_file() {
  FILE *fptr;
  fptr = fopen(OUT_FILE, "w");
  fprintf(fptr, "temp,<m>,<m2>,<E>,<E2>\n");
  fclose(fptr);
}

void write_stats(MeanStats *ms, double temp) {
  FILE *fptr;
  fptr = fopen(OUT_FILE, "a");
  fprintf(fptr, "%f,%f,%f,%f,%f\n", temp, ms->m, ms->m2, ms->E, ms->E2);
  fclose(fptr);
}

void temp_loop_simulation(Grid2D *pg, CurStats *cs) {
  int direction = END_TEMP - START_TEMP;
  MeanStats ms;

  if (direction > 0) {
    for (float temp = START_TEMP; temp < END_TEMP; temp += TEMP_STEP) {
      init_mean_stats(&ms, pg);
      simulate_temp(pg, &ms, cs, temp);
      write_stats(&ms, temp);
    }
  } else {
    for (float temp = START_TEMP; temp > END_TEMP; temp -= TEMP_STEP) {
      init_mean_stats(&ms, pg);
      simulate_temp(pg, &ms, cs, temp);
      write_stats(&ms, temp);
    }
  }
}

int main() {
  // Init rng
  srand(time(NULL)); //   seeds rand func

  // Init grid
  Grid2D grid;
  init_grid(&grid, GRID_SIZE, EXTERNAL_FIELD, INTERACTION_STRENGTH);

  // Init Stats
  CurStats cs;
  init_cur_stats(&cs, &grid);
  // Init File
  init_file();
  // Simulate
  temp_loop_simulation(&grid, &cs);

  return 0;
}
