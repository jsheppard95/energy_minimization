#Exercise 2 template for CHE210D


#import python modules
import matplotlib.pyplot as plt
import numpy as np

#import compiled fortran library for this exercise
import ex2lib

#NOTE:
#Everything below assumes unit atomic masses,
#such that forces = accelerations.


def LineSearch(Pos, Dir, dx, EFracTol,
               Accel = 1.5, MaxInc = 10., MaxIter = 10000):
    """Performs a line search along direction Dir.
Input:
    Pos: starting positions, (N,3) array
    Dir: (N,3) array of gradient direction
    dx: initial step amount, a float
    EFracTol: fractional energy tolerance
    Accel: acceleration factor
    MaxInc: the maximum increase in energy for bracketing
    MaxIter: maximum number of iteration steps
Output:
    PEnergy: value of potential energy at minimum along Dir
    Pos: minimum energy (N,3) position array along Dir
"""
    #start the iteration counter
    Iter = 0

    #find the normalized direction
    NormDir = np.clip(Dir, -1.e100, 1.e100)
    NormDir = NormDir / np.sqrt(np.sum(NormDir * NormDir))

    #take the first two steps and compute energies    
    Dists = [0., dx]
    PEs = [ex2lib.calcenergy(Pos + NormDir * x) for x in Dists]
    
    #if the second point is not downhill in energy, back
    #off and take a shorter step until we find one
    while PEs[1] > PEs[0]:
        Iter += 1
        dx = dx * 0.5
        Dists[1] = dx
        PEs[1] = ex2lib.calcenergy(Pos + NormDir * dx)
        
    #find a third point
    Dists = Dists + [2. * dx]
    PEs = PEs + [ex2lib.calcenergy(Pos + NormDir * 2. * dx)]
    
    #keep stepping forward until the third point is higher
    #in energy; then we have bracketed a minimum
    while PEs[2] < PEs[1]:
        Iter += 1
            
        #find a fourth point and evaluate energy
        Dists = Dists + [Dists[-1] + dx]
        PEs = PEs + [ex2lib.calcenergy(Pos + NormDir * Dists[-1])]

        #check if we increased too much in energy; if so, back off
        if (PEs[3] - PEs[0]) > MaxInc * (PEs[0] - PEs[2]):
            PEs = PEs[:3]
            Dists = Dists[:3]
            dx = dx * 0.5
        else:
            #we found a good energy; shift all of the points over
            PEs = PEs[-3:]
            Dists = Dists[-3:]
            dx = dx * Accel
            
    #we've bracketed a minimum; now we want to find it to high accuracy
    OldPE3 = 1.e300
    #loop over successive narrowing of the distance range
    while True:
        Iter += 1
        if Iter > MaxIter:
            print("Warning: maximum number of iterations reached in line "
                  "search.")
            break
            
        #unpack distances for ease of code-reading
        d0, d1, d2 = Dists
        PE0, PE1, PE2 = PEs

        #use a parobolic approximation to estimate the minimum location
        d10 = d0 - d1
        d12 = d2 - d1
        Num = d12*d12*(PE0-PE1) - d10*d10*(PE2-PE1)
        Dem = d12*(PE0-PE1) - d10*(PE2-PE1)
        if Dem == 0:
            #parabolic extrapolation won't work; set new dist = 0 
            d3 = 0
        else:
            #location of parabolic minimum
            d3 = d1 + 0.5 * Num / Dem
            
        #compute the new potential energy
        PE3 = ex2lib.calcenergy(Pos + NormDir * d3)
        
        #sometimes the parabolic approximation can fail;
        #check if d3 is out of range < d0 or > d2 or the new energy is higher
        if d3 < d0 or d3 > d2 or PE3 > PE0 or PE3 > PE1 or PE3 > PE2:
            #instead, just compute the new distance by bisecting two
            #of the existing points along the line search
            if abs(d2 - d1) > abs(d0 - d1):
                d3 = 0.5 * (d2 + d1)
            else:
                d3 = 0.5 * (d0 + d1)
            PE3 = ex2lib.calcenergy(Pos + NormDir * d3)
            
        #decide which three points to keep; we want to keep
        #the three that are closest to the minimum
        if d3 < d1:
            if PE3 < PE1:
                #get rid of point 2
                Dists, PEs = [d0, d3, d1], [PE0, PE3, PE1]
            else:
                #get rid of point 0
                Dists, PEs = [d3, d1, d2], [PE3, PE1, PE2]
        else:
            if PE3 < PE1:
                #get rid of point 0
                Dists, PEs = [d1, d3, d2], [PE1, PE3, PE2]
            else:
                #get rid of point 2
                Dists, PEs = [d0, d1, d3], [PE0, PE1, PE3]
                
        #check how much we've changed
        if abs(OldPE3 - PE3) < EFracTol * abs(PE3):
            #the fractional change is less than the tolerance,
            #so we are done and can exit the loop
            break
        OldPE3 = PE3

    #return the position array at the minimum (point 1)        
    PosMin = Pos + NormDir * Dists[1]
    PEMin = PEs[1]
            
    return PEMin, PosMin

        
def ConjugateGradient(Pos, dx, EFracTolLS, EFracTolCG):
    """Performs a conjugate gradient search.
Input:
    Pos: starting positions, (N,3) array
    dx: initial step amount
    EFracTolLS: fractional energy tolerance for line search
    EFracTolCG: fractional energy tolerance for conjugate gradient
Output:
    PEnergy: value of potential energy at minimum
    Pos: minimum energy (N,3) position array 
"""
    # Perform initial search in direction of forces
    forces_i_1 = ex2lib.calcforces(Pos)
    dir_i_1 = np.copy(forces_i_1)
    _, PosMin_i_1 = LineSearch(Pos, forces_i_1, dx, EFracTolLS)
    # initialize the stopping condition (dE/E) to something large
    EFracCG = 10
    while EFracCG > EFracTolCG:
        # Get energy forces at the current (last) position
        energy_i_1, forces_i = ex2lib.calcenergyforces(PosMin_i_1)
        # Set the next search direction
        # d_i = f_i + gamma_i*d_{i-1}
        gamma_i = np.dot(forces_i.flatten(), forces_i.flatten()) / \
            np.dot(forces_i_1.flatten(), forces_i_1.flatten())
        dir_i = forces_i + gamma_i*dir_i_1
        PEMin_i, PosMin_i = LineSearch(PosMin_i_1, forces_i, dx, EFracTolLS)
        EFracCG = np.abs(energy_i_1 - PEMin_i)/np.abs(PEMin_i)
        PosMin_i_1 = PosMin_i
        dir_i_1 = dir_i

    # compute final PE and return
    PE = ex2lib.calcenergy(PosMin_i)
    Pos = PosMin_i

    return PE, Pos


def InitPositions(N, L):
    """Returns an array of initial positions of each atom,
placed randomly within a box of dimensions L.
Input:
    N: number of atoms
    L: box width
Output:
    Pos: (N,3) array of positions
""" 
    Pos = np.random.uniform(high=L, size=(N, 3))
    return Pos


def KConjugateGradient(K):
    # Looping ConjugateGradient of different N's
    Nmin = 2  # minimum number of particles
    Nmax = 26  # maximum number of particles
    N_vals = np.arange(Nmin, Nmax + 1)
    # Energy Stopping conditions
    E_FRAC_TOL_LS = 1e-8
    E_FRAC_TOL_CG = 1e-10
    DX = 0.001  # LineSearch step size
    DENSITY = 0.001
    PEMins = np.zeros((len(N_vals), K))
    for i in range(len(N_vals)):
        print("N:", N_vals[i])
        # Find L s.t N/V = DENSITY
        BoxL = (N_vals[i]/DENSITY)**(1/3)
        # Perform K minimizations, each starting from different initial coords
        for j in range(K):
            # Initialize particle coordinates
            InitPos = InitPositions(N_vals[i], BoxL)
            # Get minimium energy using conjugate gradient
            PE = ConjugateGradient(InitPos, DX, E_FRAC_TOL_LS,
                                   E_FRAC_TOL_CG)[0]
            PEMins[i, j] = PE
            
    # Display results and write to file: min, average, and max of K trials for
    # each N
    min_PEs = PEMins.min(axis=1)
    mean_PEs = PEMins.mean(axis=1)
    max_PEs = PEMins.max(axis=1)
    outfile_name = f"K{K}_energy_min.txt"
    with open(outfile_name, "w") as f:
        for i in range(len(N_vals)):
            output_str = (f"N: {N_vals[i]}, "
                          f"Minimum P.E: {min_PEs[i]}, "
                          f"Average P.E: {mean_PEs[i]}, "
                          f"Maximum P.E: {max_PEs[i]}")
            print(output_str)
            f.write(output_str + "\n")


if __name__ == "__main__":
    K_vals = [100, 1000, 10000]
    for K in K_vals:
        print("K:", K)
        KConjugateGradient(K)
