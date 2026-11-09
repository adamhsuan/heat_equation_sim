import numpy as np
dt = 0.1
dl = 0.1
#L is the total spacial length of the simulation (LxLxL cube)
L = 5
K = 0.01
TOTAL_SIM_TIME = 100
NUM_LENGTH_POINTS = int(L/dl)
NUM_SIM_STEPS = int(TOTAL_SIM_TIME/dt)
REC_PRISM_TEMP = 50.0
BACKGROUND_TEMP = 20.0
REC_PRISM_X_LEN = 2
REC_PRISM_Y_LEN = 3
REC_PRISM_Z_LEN = 1

def point_in_rec_prism(point):
    if np.abs(point[0]) <= REC_PRISM_X_LEN/2 and np.abs(point[1]) <= REC_PRISM_Y_LEN/2 and np.abs(point[2]) <= REC_PRISM_Z_LEN/2:
        return True
    return False

def create_rec_prism_vtk():
    with open("rec_prism.vtk","w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("rec prism\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS 2 2 2\n")
        f.write(f"ORIGIN {-1*REC_PRISM_X_LEN/2} {-1*REC_PRISM_Y_LEN/2} {-1*REC_PRISM_Z_LEN/2}\n")
        f.write(f"SPACING {REC_PRISM_X_LEN} {REC_PRISM_Y_LEN} {REC_PRISM_Z_LEN}\n")
        
def create_temp_points_vtk(temps):

    for sim_step in range(0, NUM_SIM_STEPS):
        print(f"step {sim_step} of {NUM_SIM_STEPS}")
        with open(f"time_series_data/temperature_points_2_{sim_step}.vtk","w") as f:
            f.write('# vtk DataFile Version 3.0\n')
            f.write('temperature points\n')
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_POINTS\n")
            f.write(f"DIMENSIONS {NUM_LENGTH_POINTS+1} {NUM_LENGTH_POINTS+1} {NUM_LENGTH_POINTS+1}\n")
            f.write(f"ORIGIN {-1*L/2} {-1*L/2} {-1*L/2}\n")
            f.write(f"SPACING {dl} {dl} {dl}\n")

            f.write(f"CELL_DATA {NUM_LENGTH_POINTS**3}\n")
            temps = get_next_temps(temps)

            f.write(f"SCALARS temperature float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for i in range(NUM_LENGTH_POINTS):
                for j in range(NUM_LENGTH_POINTS):
                    for k in range(NUM_LENGTH_POINTS):
                        f.write(f"{temps[k][j][i]}\n")

#use eulers method to get the temperatures at the next time step
def get_next_temps(temps):
    next_temps = temps
    for i in range(1, NUM_LENGTH_POINTS-1):
        for j in range(1, NUM_LENGTH_POINTS-1):
            for k in range(1, NUM_LENGTH_POINTS-1):
                if point_in_rec_prism(np.array([i*dl-L/2, j*dl-L/2, k*dl-L/2])):
                    next_temps[i][j][k] = REC_PRISM_TEMP
                else:
                    next_temps[i][j][k] += K*(temps[i+1][j][k] - 2 * temps[i][j][k] + temps[i-1][j][k] + temps[i][j+1][k] - 2 * temps[i][j][k] + temps[i][j-1][k] + temps[i][j][k+1] - 2 * temps[i][j][k] + temps[i][j][k-1])/(dl**2)*dt
    return next_temps

def initiate_temps():
    initial_temps = []
    for i in range(NUM_LENGTH_POINTS):
        initial_temps.append([])
        for j in range(NUM_LENGTH_POINTS):
            initial_temps[i].append([])
            for k in range(NUM_LENGTH_POINTS):
                initial_temps[i][j].append(BACKGROUND_TEMP)
    return initial_temps

create_rec_prism_vtk()
initial_temps = np.full((NUM_LENGTH_POINTS, NUM_LENGTH_POINTS, NUM_LENGTH_POINTS), BACKGROUND_TEMP)
create_temp_points_vtk(initial_temps)
