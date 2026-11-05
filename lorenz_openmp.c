#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

// Set number of threads
#define NUM_THREADS 10

// Standard Lorenz System constants
#define SIGMA 10.0
#define RHO 28.0
#define BETA (8.0 / 3.0)

// Simulation Parameters
#define NUM_TRAJECTORIES 50 // Num of independent simulations (parallel units)
#define NUM_STEPS 10000000  // Num of time steps per trajectory
#define DT 0.00001          // Time step size
#define CHAOS 0.0001        // Value to introduce chaos

// Helps control density of information of graph
#define SAVE_EVERY 10

// Derived values
#define POINTS_PER_TRAJECTORY (NUM_STEPS / SAVE_EVERY)
#define TOTAL_POINTS (NUM_TRAJECTORIES * POINTS_PER_TRAJECTORY)

// Struct to store 3d point
typedef struct {
    double x;
    double y;
    double z;
} Point;

// Compute derivatives from Lorenz system
void lorenz_derivative(double x, double y, double z, double *dx, double *dy, double *dz) {
    *dx = SIGMA * (y - x);
    *dy = x * (RHO - z) - y;
    *dz = x * y - BETA * z;
}

// Execute one RK4 integration step
void rk4_step(double *x, double *y, double *z) {
    // k_1, k_2, k_3, k_4 slope estimates
    double k1x, k1y, k1z;
    double k2x, k2y, k2z;
    double k3x, k3y, k3z;
    double k4x, k4y, k4z;

    // Temp variables for intermediate calculations
    double temp_x, temp_y, temp_z;

    // k1: slope at current point
    lorenz_derivative(*x, *y, *z, &k1x, &k1y, &k1z);

    // k2: slope at midpoint using k1
    temp_x = *x + 0.5 * DT * k1x;
    temp_y = *y + 0.5 * DT * k1y;
    temp_z = *z + 0.5 * DT * k1z;
    lorenz_derivative(temp_x, temp_y, temp_z, &k2x, &k2y, &k2z);

    // k3: additional midpoint estimate using k2
    temp_x = *x + 0.5 * DT * k2x;
    temp_y = *y + 0.5 * DT * k2y;
    temp_z = *z + 0.5 * DT * k2z;
    lorenz_derivative(temp_x, temp_y, temp_z, &k3x, &k3y, &k3z);

    // k4: slope at end using k3
    temp_x = *x + DT * k3x;
    temp_y = *y + DT * k3y;
    temp_z = *z + DT * k3z;
    lorenz_derivative(temp_x, temp_y, temp_z, &k4x, &k4y, &k4z);

    // Combine slopes to update position
    *x += (DT / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    *y += (DT / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    *z += (DT / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);
}

// Simulate Lorenz Trajectory
void simulate_trajectory(int trajectory_id, Point *points) {
    double x = 1.0 + CHAOS * trajectory_id; // CHAOS is our value to show the chaos over time
    double y = 1.0;
    double z = 1.0;

    // Track points for current trajectory
    int save_index = 0;

    // Time evolution loop
    for (int step = 0; step < NUM_STEPS; step++) {
        // Advance system using RK4
        rk4_step(&x, &y, &z);

        // Save point based on desired save timing
        if (step % SAVE_EVERY == 0) {
            // Compute global index for current trajectory
            int index = trajectory_id * POINTS_PER_TRAJECTORY + save_index;

            // Store current state in array
            points[index].x = x;
            points[index].y = y;
            points[index].z = z;

            // Onto the next save slot for current trajectory
            save_index++;
        }
    }
}

// Write collected points to CSV file
void write_csv(const char *filename, Point *points) {
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        printf("Error: could not open output file.\n");
        return;
    }

    // CSV header
    fprintf(file, "x,y,z\n");

    // Write points
    for (int i = 0; i < TOTAL_POINTS; i++) {
        fprintf(file, "%lf,%lf,%lf\n", points[i].x, points[i].y, points[i].z);
    }

    fclose(file);
}

// Driver
int main() {
    // Allocate memory for trajectory points
    Point *points = malloc(sizeof(Point) * TOTAL_POINTS);

    if (points == NULL) {
        printf("Error: memory allocation failed.\n");
        return 1;
    }

    // Display configuration info
    printf("Lorenz Equation OpenMP Simulation\n");
    printf("Trajectories: %d\n", NUM_TRAJECTORIES);
    printf("DT: %f\n", DT);
    printf("Steps per trajectory: %d\n", NUM_STEPS);
    printf("Total saved points: %d\n", TOTAL_POINTS);
    printf("Max OpenMP threads available: %d\n", omp_get_max_threads());
    printf("Number of threads using: %d\n", NUM_THREADS);

    // Start timing
    double start_time = omp_get_wtime();

    // Adjust number of threads
    omp_set_num_threads(NUM_THREADS);

    // Parallel loop: each iteration runs num_trajectories
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < NUM_TRAJECTORIES; i++) {
        simulate_trajectory(i, points);
    }

    // Stop timing and calculate elapsed time
    double end_time = omp_get_wtime();
    double elapsed = end_time - start_time;
    printf("Simulation time: %f seconds\n", elapsed);

    // Write results and display csv file name
    write_csv("lorenz_points.csv", points);
    printf("Output written to lorenz_points.csv\n");

    // Free memory
    free(points);

    return 0;
}
